import os
import numpy as np
import pandas as pd
from intervaltree import IntervalTree, Interval
from pysam import TabixFile, AlignmentFile, FastaFile
from tqdm import tqdm
from .short_read import Coverage
from ._utils import junctions_from_cigar,splice_identical,is_same_gene, overlap
from .gene import Gene
import copy
import logging
logger=logging.getLogger('isotools')

#### io functions for the transcriptome class

def add_short_read_coverage(self, bam_files,names=None, load=False):
    'Add short read coverage to the genes.\n This does, by default (e.g. if load==False), not actually read the bams, but reading is done at first access'
    self.infos.setdefault('short_reads', pd.DataFrame(columns=['name', 'file']))
    if isinstance(bam_files,str) and isinstance(names,str):
        bam_files=[bam_files]
        names=[names]
    if isinstance(bam_files,list) and isinstance(names, list):
        assert len(bam_files) ==len(names), 'lists of sample names and bam files must have same length'
    elif isinstance(bam_files, dict):
        names=list(bam_files)
        bam_files=list(bam_files.values())
    else:
        raise ValueError('either provide a bam file and names as strings, lists of strings or as a dict')
    self.infos['short_reads']=pd.concat([self.infos['short_reads'], pd.DataFrame({'name':names, 'file':bam_files})])
    if load: # when loading coverage for all genes keep the filehandle open, hopefully a bit faster
        for i,bamfile in enumerate(self.infos['short_reads'].file):            
            logger.info('Adding short read coverag from %s',bamfile)
            if i==0:
                for g in tqdm(self):
                    g.data.setdefault('short_reads',list())
            with AlignmentFile(bamfile, "rb") as align:
                for g in tqdm(self):
                    if len(g.data['short_reads'])==i:
                        g.data['short_reads'].append(Coverage.from_alignment(align,g))
                
 

def remove_samples(self, sample_names):
    ''' removes samples from the dataset'''
    if isinstance(sample_names, str):
        sample_names=[sample_names]
    assert all(s in self.samples for s in sample_names), 'Did not find all samples to remvoe in dataset'
    rm_idx=self.sample_table.index[self.sample_table.name.isin(sample_names)]
    self.sample_table.drop(index=rm_idx)
    for g in self:
        remove_tr=[]
        for i,tr in enumerate(g.transcripts):
            if any(s in tr['samples'] for s in sample_names):
                tr['samples']={s:d for s,d in tr['samples'].items() if s not in sample_names}
                if not tr['samples']:
                    remove_tr.append(i)
        if remove_tr: #remove the transcripts that is not expressed by remaining samples
            g.data['transcripts']=[tr for i,tr in enumerate(g.transcripts) if i not in remove_tr]
            g.data['splice_graph']=None # gets recomputed on next request
            
        elif 'splice_graph' in g.data:
            g.splice_graph.weights=np.delete(g.splice_graph.weights,rm_idx, 0)

def add_sample_from_bam(self,fn, sample_name,fuzzy_junction=5,add_chromosomes=True,  **kwargs):
    '''import expressed transcripts from bam and add it to existing transcriptome'''
    #todo: one alignment may contain several samples - this is not supported at the moment
    assert sample_name not in self.samples, 'sample %s is already in the data set.' % sample_name
    logger.info(f'adding sample {sample_name}')
    kwargs['name']=sample_name
    kwargs['file']=fn
    #genome_fh=FastaFile(genome_fn) if genome_fn is not None else None
    with AlignmentFile(fn, "rb") as align:        
        if add_chromosomes:
            chromosomes=align.references
        else:
            chromosomes=self.chromosomes
        stats = align.get_index_statistics()
        # try catch if sam/ no index /not pacbio?
        total_alignments = sum([s.mapped for s in stats if s.contig in chromosomes])
        total_reads=0
        chimeric=dict()
        with tqdm(total=total_alignments, unit='transcripts') as pbar:
            for chrom in chromosomes: #todo: potential issue here - secondary/chimeric alignments to non listed chromosomes are ignored
                pbar.set_postfix(chr=chrom)
                #transcripts=IntervalTree()
                #novel=IntervalTree()
                chr_len=align.get_reference_length(chrom)
                transcripts=IntervalArray(chr_len)
                novel=IntervalArray(chr_len)
                for read in align.fetch(chrom):
                    pbar.update()
                    if read.flag & 0x100:
                        n_secondary+=1
                        continue # use only primary alignments
                    strand = '-' if read.is_reverse else '+'
                    exons = junctions_from_cigar(read.cigartuples,read.reference_start)
                    tags= dict(read.tags)
                    tr_range=(exons[0][0], exons[-1][1])
                    if 'is' in tags:
                        cov=tags['is'] #number of actual reads supporting this transcript
                    else:
                        cov=1
                    total_reads+=cov
                    for tr_interval in transcripts.overlap(*tr_range):
                        tr=tr_interval.data
                        if tr['strand'] != strand:
                            continue
                        if splice_identical(exons, tr['exons']):                                
                            tr['samples'].setdefault(sample_name,{}).setdefault('range',{}).setdefault(tr_range,0)
                            tr['samples'][sample_name]['range'][tr_range]+=cov                               
                            break
                    else:
                        tr={'exons':exons,'samples':{sample_name:{'range':{tr_range:cov}}}, 'strand':strand}
                        transcripts.add(Interval(*tr_range,tr))
                    #if genome_fh is not None:
                    #    mutations=get_mutations(read.cigartuples, read.query_sequence, genome_fh, chrom,read.reference_start,read.query_qualities)
                    #    for pos,ref,alt,qual in mutations:
                    #        tr['samples'][sample_name].setdefault('mutations',{}).setdefault(pos,{'ref':ref}).setdefault(alt,[0,[]])
                    #        tr['samples'][sample_name]['mutations'][pos][alt][0]+=cov
                    #        if qual:
                    #            tr['samples'][sample_name]['mutations'][pos][alt][1].append(qual) #assuming the quality already accounts for cov>1
                    if 'SA' in tags:#part of a chimeric alignment
                        chimeric.setdefault(read.query_name, []).append({'chr':chrom,'strand':strand,'range':tr_range,'cigar':read.cigarstring, 'aligned': aligned_part(read.cigartuples, read.is_reverse)}) #cigar is required to sort out the order 
                        tr['samples'][sample_name].setdefault('chimeric',[]).append(chimeric[read.query_name])                       
                    elif 4 in read.cigartuples: #clipping
                        clip=get_clipping(read.cigartuples, read, read.reference_start)
                        tr['samples'][sample_name].setdefault('clipping',{}).setdefault(clip,0)
                        tr['samples'][sample_name]['clipping'][clip]+=cov
                for tr_interval in transcripts:
                    tr=tr_interval.data
                    tr_ranges=tr['samples'][sample_name].pop('range')
                    tr['exons'][0][0]=min(r[0] for r in tr_ranges) #todo - instead of extreams take the median?
                    tr['exons'][-1][1]=max(r[1] for r in tr_ranges)  
                    tr['samples'][sample_name].setdefault('coverage',0)
                    tr['samples'][sample_name]['coverage']+=sum(tr_ranges.values())                       
                    gene=self._add_sample_transcript(tr,chrom, sample_name, fuzzy_junction) #tr is not updated
                    if gene is None:
                        novel.add(tr_interval)    
                self._add_novel_genes(novel,chrom)
    kwargs['total_reads']=total_reads
    self.infos['sample_table']=self.sample_table.append(kwargs, ignore_index=True)

    for g in self:
        if 'splice_graph' in g.data and g.data['splice_graph'] is not None: # still valid splice graphs no new transcripts - add a row of zeros to weights
            n_samples=len(self.samples)
            w=g.coverage
            if w.shape[0]<n_samples:                    
                w_new=np.zeros((n_samples,w.shape[1]))
                w_new[:w.shape[0],:] = w
                g.data['splice_graph'].weights=w_new
    self.make_index()#assure the index is correct

def _add_sample_transcript(self,tr,chrom,sample_name,fuzzy_junction=5):
    'add transcript to gene in chrom - return gene on success and None if no Gene was found'
    if chrom not in self.data:
        return None
    genes_ol = [g for g in self.data[chrom][tr['exons'][0][0]: tr['exons'][-1][1]] if g.strand==tr['strand'] ] 
    #check if transcript is already there (e.g. from other sample):
    for g in genes_ol:
        for tr2 in g.transcripts:
            if splice_identical(tr2['exons'], tr['exons']):
                tr2['samples'][sample_name]=tr['samples'][sample_name]
                return g
    #check if gene is already there (e.g. from same or other sample):
    g,(sj_i, base_i)=_get_intersects(genes_ol, tr['exons'])
    if g is not None:            
        if g.is_annotated: #check for fuzzy junctions (e.g. same small shift at 5' and 3' compared to reference annotation)
            shifts=g.correct_fuzzy_junctions(tr, fuzzy_junction, modify=True) #this modifies tr['exons']
            if shifts: 
                tr['samples'][sample_name]['fuzzy_junction']=shifts #keep the info, mainly for testing/statistics    
                for tr2 in g.transcripts:#check if correction made it identical to existing
                    if splice_identical(tr2['exons'], tr['exons']):
                        tr2['samples'][sample_name]=tr['samples'][sample_name]
                        return g 
            altsplice= g.ref_splice_graph.get_alternative_splicing(tr['exons'], g.strand)
            tr['annotation']={'sj_i': sj_i, 'base_i':base_i,'as':altsplice}

        else: # range of the novel gene might have changed
            start, end=min(tr['exons'][0][0],g.start),max(tr['exons'][-1][1],g.end)
            if start<g.start or end>g.end: #todo: potential issue: in this case two genes may have grown together
                new_gene=Gene(start, end,g.data,self)
                self.data[chrom].add(new_gene)
                self.data[chrom].remove(g)
                g=new_gene
        g.data.setdefault('transcripts',[]).append(tr) 
        g.data['splice_graph']=None #gets recomputed on next request           
    return g

def _add_novel_genes( self,novel,chrom, spj_iou_th=0, reg_iou_th=.5, gene_prefix='PB_novel_'):
    '"novel" is a tree of transcript intervals (not Gene objects) ,e.g. from one chromosome, that do not overlap any annotated or unanntoated gene'
    n_novel=self.novel_genes
    idx={id(tr):i for i,tr in enumerate(novel)}
    merge=list()
    for i,tr in enumerate(novel):    
        merge.append({tr})
        candidates=[c for c in novel.overlap(tr.begin, tr.end) if c.data['strand']==tr.data['strand'] and idx[id(c)]<i]
        for c in candidates:
            if c in merge[i]:
                continue
            if is_same_gene(tr.data['exons'], c.data['exons'],spj_iou_th, reg_iou_th):
                #add all transcripts of candidate
                merge[i].update(merge[idx[id(c)]])
        for c in merge[i]:#update all overlapping (add the current to them)
            merge[idx[id(c)]]=merge[i]
            
    seen=set()
    for trS in merge: 
        if id(trS) in seen:
            continue
        seen.add(id(trS))
        trL=[tr.data for tr in trS]
        strand=trL[0]['strand']
        start=min(tr['exons'][0][0] for tr in trL)
        end=max(tr['exons'][-1][1] for tr in trL)
        n_novel+=1
        new_data={'chr':chrom, 'ID':f'{gene_prefix}{n_novel:05d}', 'strand':strand, 'transcripts':trL}
        self.data.setdefault(chrom,IntervalTree()).add(Gene(start,end,new_data,self ))
    self.infos['novel_counter']=n_novel

def _get_intersects(genes_ol, exons):
    'calculate the intersects of all overlapping genes and return the best'  
    ## todo: in order to detect readthrough fusion genes, report all overlapping reference genes? or exon/junction wise?
    # prefer annotated genes
    intersect=[(i,g.ref_splice_graph.get_intersects(exons)) for i,g in enumerate(genes_ol) if g.is_annotated]
    if intersect:
        best_idx ,intersects= max(intersect, key=lambda x:x[1])
        if intersects[1] > 0 :    # at least some overlap (definition could be changed here)
            return genes_ol[best_idx],intersects
    #no annotated overlap, look for novel genes
    intersect=[(i,g.splice_graph.get_intersects(exons)) for i,g in enumerate(genes_ol) if not g.is_annotated]
    if intersect:
        best_idx ,intersects= max(intersect, key=lambda x:x[1])
        if intersects[1] > 0 :    # at least some overlap (definition could be changed here)
            return genes_ol[best_idx],intersects
    #todo: potentially antisense could be considered as well here
    return None, (0,0)

def import_gtf_transcripts(fn,transcriptome, chromosomes=None):
    '''import transcripts from gtf file (e.g. for a reference)
    returns a dict interval trees for the genes, each containing the splice graph'''
    gtf = TabixFile(fn)   
    exons = dict()  # transcript id -> exons
    transcripts = dict()  # gene_id -> transcripts
    skipped = set()
    genes=dict()  
    gene_dict={}
    for line in tqdm(gtf.fetch(), smoothing=.1): 
        ls = line.split(sep="\t")        
        if chromosomes is not None and ls[0] not in chromosomes:
            #warnings.warn('skipping line from chr '+ls[0])
            continue
        _=genes.setdefault(ls[0], IntervalTree())
        info = dict([pair.lstrip().split(' ', 1) for pair in ls[8].replace('"','').split(";") if pair])        
        start, end = [int(i) for i in ls[3:5]]
        start -= 1  # to make 0 based
        if ls[2] == "exon":
            # logger.debug(line)
            try:
                gtf_id = info['transcript_id']
                _=exons.setdefault(gtf_id, list()).append((start, end))
            except KeyError:  # should not happen if GTF is OK
                logger.warn("gtf format error: exon without transcript_id. Skipping line\n"+line)
        elif ls[2] == 'gene':
            if 'gene_id' not in info:
                logger.warn("gtf format error: gene without gene_id. Skipping line\n"+line)
            else:
                info=_prepare_gene_info(info, ls[0], ls[6])
                new_gene=Gene(start, end, info, transcriptome)
                genes[ls[0]].add(new_gene)
                gene_dict[info['gene_id']]=new_gene
        elif ls[2] == 'transcript':
            try:
                _=transcripts.setdefault(info["gene_id"], list()).append(info["transcript_id"])
            except KeyError:
                logger.warn("gtf format errror: transcript must have gene_id and transcript_id, skipping line\n"+line )
      
        else:
            skipped.add(ls[2])
    if len(skipped)>0:
        logger.info('skipped the following categories: {}'.format(skipped))
    # sort the exons
    for tid in exons.keys():
        exons[tid].sort()    
    #check for missing gene information
    missed_genes=0
    for gid in transcripts:
        if gid not in gene_dict:
            missed_genes+=1 #alternativly add a gene with name of transcript
    if missed_genes:
        raise ValueError('no specific gene information for {}/{} genes'.format(missed_genes,missed_genes+sum((len(t) for t in genes.values()) )))
    
    for chrom in genes:#link transcripts to genes
        for gene in genes[chrom]:
            g_id = gene.id
            t_ids = transcripts.get(g_id,  g_id)
            gene.data['reference'].setdefault('transcripts', [])
            for t_id in t_ids:
                try:
                    gene.data['reference']['transcripts'] = [{'transcript_id': t_id, 'exons':exons[t_id]}]
                except KeyError:
                    # genes without transcripts get a single exons transcript
                    gene.data['reference']['transcripts'] =  [{'transcript_id': t_id,'exons':[tuple(gene[:2])]}]
    return genes,{'reference_file':fn}

def import_gff_transcripts(fn, transcriptome, chromosomes=None, gene_categories=['gene']):
    '''import transcripts from gff file (e.g. for a reference)
    returns a dict interval trees for the genes'''
    file_size=os.path.getsize(fn)
    gff = TabixFile(fn)        
    chrom_ids = get_gff_chrom_dict(gff, chromosomes)    
    exons = dict()  # transcript id -> exons
    transcripts = dict()  # gene_id -> transcripts
    skipped = set()    
    genes=dict()
    gene_set=set()
    cds_start=dict()
    cds_stop=dict()
    #takes quite some time... add a progress bar?
    for line in tqdm(gff.fetch(), smoothing=.1):  
        ls = line.split(sep="\t")
        if ls[0] not in chrom_ids:
            continue
        chrom = chrom_ids[ls[0]]
        genes.setdefault(chrom, IntervalTree())
        try:
            info = dict([pair.split('=', 1) for pair in ls[8].split(";")])
        except ValueError:
            logger.warning("GFF format error in infos (should be ; seperated key=value pairs). Skipping line: "+line)
        start, end = [int(i) for i in ls[3:5]]
        start -= 1  # to make 0 based
        if ls[2] == "exon":
            try:
                gtf_id = info['Parent']
                exons.setdefault(gtf_id, list()).append((start, end))
            except KeyError:  # should not happen if GTF is OK
                logger.warning("GFF format error: no parent found for exon. Skipping line: "+line)
        elif ls[2] == 'gene' or 'ID' in info and info['ID'].startswith('gene'):
            info['strand'] = ls[6]
            info['chr'] = chrom
            #genes[chrom][start:end] = info
            _set_alias(info,{'ID':['gene_id'],'name':['Name','gene_name']})
            gene_set.add(info['ID'])
            
            ref_info={k:v for k,v in info.items() if k not in Gene.required_infos}
            info={k:info[k] for k in Gene.required_infos}
            info['reference']=ref_info
            genes[chrom].add(Gene(start,end,info, transcriptome))
        elif all([v in info for v in ['Parent', "ID"]]) and (ls[2] == 'transcript' or info['Parent'].startswith('gene')):# those denote transcripts
            tr_info={k:v for k,v in info.items() if k.startswith('transcript_')}
            transcripts.setdefault(info["Parent"], {})[info["ID"]]=tr_info
        elif ls[2]  == 'start_codon' and 'Parent' in info:
            cds_start[info['Parent']]=end if ls[6]=='-' else start
        elif ls[2]  == 'stop_codon' and 'Parent' in info:
            cds_stop[info['Parent']]=start if ls[6]=='-' else end            
        else:
            skipped.add(ls[2]) #transcript infos?
    if skipped:
        logger.info('skipped the following categories: {}'.format(skipped))
    # sort the exons
    logger.info('sorting exon positions...')
    for tid in tqdm(exons):
        exons[tid].sort()
    missed_genes=[gid for gid in transcripts.keys() if gid not in gene_set]
    if missed_genes:        
        #logger.debug('/n'.join(gid+str(tr) for gid, tr in missed_genes.items()))
        notfound=len(missed_genes)
        found=sum((len(t) for t in genes.values()) )
        logger.warning('Missing genes! Found gene information in categories {} for {}/{} genes'.format(gene_categories,found, found+notfound))
    logger.info('building gene data structure...')
    # add transcripts to genes
    for chrom in genes:
        for gene in genes[chrom]:
            
            g_id = gene.id
            tr=transcripts.get(g_id,  {g_id:{}})
            for t_id,tr_info in tr.items():
                tr_info['transcript_id']=t_id
                try:
                    tr_info['exons']=exons[t_id]
                except KeyError:
                    # genes without exons get a single exons transcript
                    tr_info['exons']=[tuple(gene[:2])]
                #add cds
                if t_id in cds_start and t_id in cds_stop:
                    tr_info['CDS']=(cds_start[t_id], cds_stop[t_id]) if cds_start[t_id]< cds_stop[t_id] else (cds_stop[t_id],cds_start[t_id])
                gene.data['reference'].setdefault('transcripts', []).append(tr_info)
    return genes,{'reference_file':fn}

## io utility functions

def get_mutations(cigartuples, seq, ref, chrom,ref_start,qual):
    'look up the bases affected by mutations as reported in the cigar string'
    #cigar numbers:
    #012345678
    #MIDNSHP=X
    mutations=[]
    ref_pos=ref_start
    seq_pos=0
    for cigar in cigartuples:
        if cigar[0] in (1,2,8):#I(ins), D(del) or X (missmatch): 
            ref_base='' if cigar[0] == 1 else ref.fetch(chrom, ref_pos,ref_pos+cigar[1])
            alt_base='' if cigar[0] == 2 else seq[seq_pos:(seq_pos+cigar[1])]
            mutations.append((ref_pos,ref_base, alt_base,qual[seq_pos] if qual else None))
        if cigar[0] in (0, 2, 3, 7, 8):  # MDN=X -> move forward on reference
            ref_pos += cigar[1]
        if cigar[0] in (0, 1, 4,7, 8):  # MIS=X -> move forward on seq
            seq_pos += cigar[1]
    return mutations      

def aligned_part(cigartuples, is_reverse):
    "returns the interval of the trasncript that is aligned (e.g. not clipped) according to cigar. Positions are according to transcript strand"
    parts=[]
    start=end=0
    for cigar in reversed(cigartuples) if is_reverse else cigartuples:
        if cigar[0] in (0, 1, 7, 8):  # MI=X -> move forward on read:
            end+=cigar[1]
        elif cigar[0] in (4,5): #clipping at end
            if end>start:
                return (start,end)
            end+=cigar[1]
            start=end
    return (start, end) #clipping at begining or no clipping

def get_clipping(cigartuples, pos, is_reverse):
    if cigartuples[0][0]==4:
        #clipping at the begining
        return(pos, -cigartuples[0][1])
    elif cigartuples[-1][0]==4:
        #clipping at the end - get the reference position
        return(pos+sum(c[1] for c in cigartuples[:-1] if c[0] in (0, 2, 3, 7, 8)), cigartuples[-1][1])# MDN=X -> move forward on reference:
    else:
        return None

def _set_alias(d,alias):
    for pref,alt in alias.items():
        alt=[a for a in alt if a in d]
        if pref not in d:
            try:
                d[pref]=next(d[a] for a in alt)
            except StopIteration:
                logger.error(f'did not find alternative for {pref}- suggested terms are {alt}, but have only those keys:{list(d)}')
                raise
        for a in alt:
            d.pop(a,None)


def _prepare_gene_info(info,chrom, strand, modify=True):
    if not modify:
        info=copy.deepcopy(info)
    _set_alias(info, {'name':['Name','gene_id']})
    assert 'ID' in info
    info['strand'] =strand
    info['chr'] = chrom
    ref_info={k:v for k,v in info.items() if k not in Gene.required_infos}
    info={k:info[k] for k in Gene.required_infos}
    info['reference']=ref_info
    return info

#human readable output
def gene_table(self, region=None ): #ideas: filter, extra_columns
    'create a gene summary table'
    colnames=['chr', 'start', 'end', 'strand', 'gene_name', 'n_transcripts']        
    rows=[(g.chrom, g.start, g.end, g.strand, g.id, g.n_transcripts) for g in  self.iter_genes(region)]
    df = pd.DataFrame(rows, columns=colnames)
    return(df)

def transcript_table(self, region=None, extra_columns=None,  include=None, remove=None): 
    'create a transcript table'
    if extra_columns is None:
        extra_columns=[]
    if not isinstance( extra_columns, list):
        raise ValueError('extra_columns should be provided as list')
    extracolnames={'annotation':('splice_intersect', 'base_intersect','splicing_comparison')}
    if 'sample_table' in self.infos:
        extracolnames['coverage']=(f'coverage_{sn}' for sn in self.infos['sample_table'].name)            
        #extracolnames['group_coverage']=(f'coverage_{sn}' for sn in self.infos['sample_table'].groups)
    colnames=['chr', 'gene_start', 'gene_end','strand', 'gene_id','gene_name' ,'transcript_name']     + [n for w in extra_columns for n in (extracolnames[w] if w in extracolnames else (w,))]
    rows=[]
    for g,trid, _ in self.iter_transcripts(region, include,remove):                
            rows.append((g.chrom, g.start, g.end, g.strand,g.id, g.name, trid)+g.get_values(trid,extra_columns))
    df = pd.DataFrame(rows, columns=colnames)
    return(df)

def fusion_table(self, region=None,  include=None, remove=None, star_chimeric=None, illu_len=200):
    'find fusion genes and compile a table with relevant infos (breakpoints, coverage, ...)'
    #todo: correct handeling of three part fusion events not yet implemented
    #todo: ambiguous alignment handling not yet implemented
    logger.warn('fusion table is broken')
    if star_chimeric is None:
        star_chimeric=dict()
    assert isinstance(star_chimeric, dict)
    
    fusion_coverage=list()
    breakpoints=dict()

    for g,trid,tr in self.iter_transcripts(region=region, include=include, remove=remove):            
        if 'fusion' in tr:
            total_cov=sum(tr['coverage'])
            greg1=(tr['exons'][0][0],tr['exons'][-1][-1])      
            pos1=greg1[int(g.strand=="+")]      
            breakpoint1=f'{g.chrom}:{pos1}({g.strand})'
            for f in tr['fusion']:
                greg2=(f['exons'][0][0],f['exons'][-1][-1])
                pos2=greg2[int(f["strand"]=="+")]
                breakpoint2=f'{f["chr"]}:{pos2}({f["strand"]})'
                gene2='intergenic' if f['annotation'] is None else f['annotation']['ref_gene_name']
                if tr['aligned'][0]<f['aligned'][0]:
                    fusion_coverage.append([trid,tr['source_len'],g.name,tr['aligned'],breakpoint1,gene2,f['aligned'],breakpoint2,total_cov]+list( tr['coverage'])+[0]*len(star_chimeric))
                else:
                    fusion_coverage.append([trid,tr['source_len'],gene2,f['aligned'],breakpoint2,g.name,tr['aligned'],breakpoint1,total_cov]+list( tr['coverage'])+[0]*len(star_chimeric))
                breakpoints.setdefault(g.chrom, IntervalTree()).addi(pos1-illu_len*int(g.strand=='+'), pos1+illu_len*int(g.strand=='-'),len(fusion_coverage)-1)
                breakpoints.setdefault(f["chr"], IntervalTree()).addi(pos2-illu_len*int(f['strand']=='+'), pos2+illu_len*int(f['strand']=='-'),len(fusion_coverage)-1)
    offset=9+len(self.infos['sample_table'])
    for sa_idx,sa in enumerate(star_chimeric):
        chimeric_tab=pd.read_csv(star_chimeric[sa], sep='\t')
        for _, row in chimeric_tab.iterrows():
            if row['chr_donorA'] in breakpoints and row['chr_acceptorB'] in breakpoints:
                idx1={bp.data for bp in breakpoints[row['chr_donorA']][row['brkpt_donorA']]}
                if idx1:
                    idx2={bp.data for bp in breakpoints[row['chr_acceptorB']][row['brkpt_acceptorB']]}
                    idx_ol={idx for idx,snd in idx2 if (idx, not snd) in idx1}
                    for idx in idx_ol:
                        fusion_coverage[idx][offset+sa_idx]+=1
                    
    fusion_coverage=pd.DataFrame(fusion_coverage, columns=['trid','len', 'gene1','part1','breakpoint1', 'gene2','part2','breakpoint2','total_cov']+[s+'_cov' for s in self.infos['sample_table'].name]+[s+"illu_cov" for s in star_chimeric])
    return fusion_coverage        

def write_gtf(self, fn, source='isotools',use_gene_name=False,  include=None, remove=None):     
        'writes the transcripts in gtf format to a file'
        with open(fn, 'w') as f:     
            for gene in tqdm(self):
                lines=gene.to_gtf(source=source,use_gene_name=use_gene_name, include=include, remove=remove)
                if lines:
                    _=f.write('\n'.join( ('\t'.join(str(field) for field in line) for line in lines) )+'\n')


def get_gff_chrom_dict(gff, chromosomes):
    'fetch chromosome ids - in case they use ids in gff for the chormosomes'
    chrom = {}
    for c in gff.contigs:
        #loggin.debug ("---"+c)
        for line in gff.fetch(c, 1, 2): #chromosomes span the entire chromosome, so they can be fetched like that
            if line[1] == "C":
                ls = line.split(sep="\t")
                if ls[2] == "region":
                    info = dict([pair.split("=")
                                    for pair in ls[8].split(";")])
                    if "chromosome" in info.keys():
                        if chromosomes is None or info["chromosome"] in chromosomes:
                            chrom[ls[0]] = info["chromosome"]
                        break
                        
        else:# no specific regions entries - no aliases
            if chromosomes is None or c in chromosomes:
                chrom[c]=c
    return(chrom)

           

class IntervalArray:
    #drop in replacement for the interval tree during construction, with faster lookup
    def __init__(self, total_size, bin_size=1e4):
        self.obj={}
        self.data=[set() for _ in range(int((total_size)//bin_size)+1)]
        self.bin_size=bin_size

    def overlap(self,begin, end):
        candidates={obj_id for idx in range(int(begin//self.bin_size), int(end//self.bin_size)+1) for obj_id in self.data[idx]}
        return (self.obj[obj_id] for obj_id in candidates if overlap((begin, end),self.obj[obj_id]) ) #this assumes object has range obj[0] to obj[1]

    def add(self,obj):
        self.obj[id(obj)]=obj
        for idx in range(int(obj.begin//self.bin_size), int(obj.end//self.bin_size)+1):
            self.data[idx].add(id(obj))
    
    def __iter__(self):
        return (v for v in self.obj.values())