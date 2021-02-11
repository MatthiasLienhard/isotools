import numpy as np
import pandas as pd
from intervaltree import IntervalTree, Interval
from pysam import TabixFile, AlignmentFile, FastaFile
from tqdm import tqdm
from .short_read import Coverage
from ._utils import junctions_from_cigar,splice_identical,is_same_gene, overlap, pairwise
from .gene import Gene
import copy
import logging
logger=logging.getLogger('isotools')

#### io functions for the transcriptome class
SPLICE_CATEGORY=['FSM','ISM','NIC','NNC','NOVEL']
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
            with AlignmentFile(bamfile, "rb") as align:
                for g in tqdm(self):
                    g.data.setdefault('short_reads',list())
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
            if any(s in tr['coverage'] for s in sample_names):
                tr['coverage']={s:cov for s,cov in tr['coverage'].items() if s not in sample_names}
                if not tr['coverage']:
                    remove_tr.append(i)
        if remove_tr: #remove the transcripts that is not expressed by remaining samples
            g.data['transcripts']=[tr for i,tr in enumerate(g.transcripts) if i not in remove_tr]
        g.data['splice_graph']=None # gets recomputed on next request
        g.data['coverage']=None            
        
def add_sample_from_bam(self,fn, sample_name,fuzzy_junction=5,add_chromosomes=True,chimeric_mincov=2,  **kwargs):
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
        with tqdm(total=total_alignments, unit='reads') as pbar:
            
            for chrom in chromosomes: #todo: potential issue here - secondary/chimeric alignments to non listed chromosomes are ignored
                pbar.set_postfix(chr=chrom)
                #transcripts=IntervalTree()
                #novel=IntervalTree()
                chr_len=align.get_reference_length(chrom)
                transcripts=IntervalArray(chr_len) #intervaltree was pretty slow for this context 
                novel=IntervalArray(chr_len)
                n_reads=0
                for read in align.fetch(chrom):
                    n_reads+=1
                    pbar.update(.5)
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
                    if 'SA' in tags:#part of a chimeric alignment
                        if chimeric_mincov>0:
                            chimeric.setdefault(read.query_name, [{sample_name:cov},[]])[1].append([chrom,strand,exons, aligned_part(read.cigartuples, read.is_reverse),None]) 
                        continue

                    for tr_interval in transcripts.overlap(*tr_range):    #did we see this transcript already?
                        if tr_interval.data['strand'] != strand:
                            continue
                        if splice_identical(exons, tr_interval.data['exons']):   
                            tr=tr_interval.data                             
                            tr.setdefault('range',{}).setdefault(tr_range,0)
                            tr['range'][tr_range]+=cov                               
                            break
                    else:
                        tr={'exons':exons,'range':{tr_range:cov}, 'strand':strand}
                        transcripts.add(Interval(*tr_range,tr))
                    #if genome_fh is not None:
                    #    mutations=get_mutations(read.cigartuples, read.query_sequence, genome_fh, chrom,read.reference_start,read.query_qualities)
                    #    for pos,ref,alt,qual in mutations:
                    #        tr.setdefault('mutations',{}).setdefault(sample_name,{}).setdefault(pos,{'ref':ref}).setdefault(alt,[0,[]])
                    #        tr['mutations'][sample_name][pos][alt][0]+=cov
                    #        if qual:
                    #            tr['mutations'][sample_name][pos][alt][1].append(qual) #assuming the quality already accounts for cov>1
                                    
                    if 4 in read.cigartuples: #clipping
                        clip=get_clipping(read.cigartuples, read, read.reference_start)
                        tr.setdefault('clipping',{}).setdefault(sample_name,{}).setdefault(clip,0)
                        tr['clipping'][sample_name][clip]+=cov
                for tr_interval in transcripts:
                    tr=tr_interval.data
                    tr_ranges=tr.pop('range')
                    tr['exons'][0][0]=min(r[0] for r in tr_ranges) #todo - instead of extreams take the median?
                    tr['exons'][-1][1]=max(r[1] for r in tr_ranges)  
                    tr.setdefault('coverage',{}).setdefault(sample_name,0)
                    cov=sum(tr_ranges.values())
                    tr['coverage'][sample_name]+=cov
                    gene=self._add_sample_transcript(tr,chrom, sample_name, fuzzy_junction) 
                    if gene is None:
                        novel.add(tr_interval)  
                    n_reads-=cov
                    pbar.update(cov/2)
                self._add_novel_genes(novel,chrom)
                pbar.update(n_reads/2) #the chimeric reads not processed here, add them to the progress
    kwargs['total_reads']=total_reads
    self.infos['sample_table']=self.sample_table.append(kwargs, ignore_index=True)
    #merge chimeric reads and assign gene names
    chimeric,non_chimeric=_check_chimeric(chimeric)
    self._add_chimeric(chimeric, chimeric_mincov)
    novel=dict()
    for cov,((chrom,strand,exons,parts, _),) in non_chimeric:
        tr={'exons':exons, 'coverage':cov,'strand':strand,'chr':chrom, 'long_intron_chimeric':[parts] }
        gene=self._add_sample_transcript(tr,chrom, sample_name, fuzzy_junction) #tr is not updated
        if gene is None:
            novel.setdefault(chrom,[]).append(tr)
    for chrom in novel:
        self._add_novel_genes(IntervalTree(Interval(tr['exons'][0][0],tr['exons'][-1][1], tr) for tr in novel[chrom]),chrom)
    #self.infos.setdefault('chimeric',{})[sample_name]=chimeric #save all chimeric reads (for debugging)
    for g in self:
        if 'coverage' in g.data and g.data['coverage'] is not None: # still valid splice graphs no new transcripts - add a row of zeros to weights
            g._set_coverage()

def _add_chimeric(self,new_chimeric, min_cov):
    ''' add new chimeric transcripts to transcriptome, if covered by > min_cov reads
    '''        
    for new_bp,new_chim_list in new_chimeric.items():        
        if sum(sum(cov.values()) for cov,_ in new_chim_list)<min_cov:
            continue
        for new_chim in new_chim_list:
            #todo: check special cases (long intron, one part only), discard invalid (large overlaps, large gaps)
            #find equivalent chimeric reads
            for found in self.chimeric.setdefault(new_bp,[]):
                if all(splice_identical(ch1[2], ch2[2]) for ch1,ch2 in zip(new_chim[1],found[1])):
                    for sa in new_chim[0]:#add coverage
                        found[0][sa]=found[0].get(sa,0)+new_chim[0][sa]
                    #adjust start
                    if found[1][0][1]=='+':#strand of first part
                        found[1][0][2][0][0]=min(found[1][0][2][0][0],new_chim[1][0][2][0][0])
                    else:
                        found[1][0][2][0][1]=max(found[1][0][2][0][1],new_chim[1][0][2][0][1])
                    #adjust end
                    if found[1][-1][1]=='+':#strand of last part
                        found[1][-1][2][-1][1]=max(found[1][-1][2][-1][1],new_chim[1][-1][2][-1][1])
                    else:
                        found[1][-1][2][-1][0]=min(found[1][-1][2][-1][0],new_chim[1][-1][2][-1][0])
                    break
            else: #not seen
                self.chimeric[new_bp].append(new_chim)
                for part in new_chim[1]:
                    if part[0] in self.data:
                        genes_ol=[g for g in self.data[part[0]][part[2][0][0]: part[2][-1][1]] if g.strand==part[1]]
                        #g,_=_get_intersects(genes_ol, part[2])
                        g,_,_,_=_find_splice_sites(genes_ol, part[2]) #take the best - ignore other hits here
                        if g is not None:
                            part[4]=g.name
                            g.data.setdefault('chimeric',{})[new_bp]=self.chimeric[new_bp]
                

def _breakpoints(chimeric):
    ''' gets chimeric aligment as a list and returns list of breakpoints.
        each breakpoint is a tuple of (chr1, strand1, pos1,  chr2,strand2,pos2)
    '''
    return tuple( (a[0],a[1],a[2][-1][1] if a[1]=='+' else a[2][0][0], 
                  b[0],b[1],b[2][0][0]  if b[1]=='+' else b[2][-1][1])
                for a,b in pairwise(chimeric))

def _check_chimeric(chimeric):
    ''' prepare the chimeric reads:
        1) sort parts according to read order
        2) compute breakpoints
        3) check if the chimeric read is actually a long introns - return list as nonchimeric
        4) sort into dict by breakpoint - return dict as chimeric

    '''
    chimeric_dict={}
    non_chimeric=list()
    for new_chim in chimeric.values():    
        new_chim[1].sort(key=lambda x: x[3][1]) #sort by end of aligned part
        bpts=_breakpoints(new_chim[1])#compute breakpoints
        merge=[i for i,bp in enumerate(bpts) if
                bp[0]==bp[3] and #same chr
                bp[1]==bp[4] and #same strand,
                0 < (bp[5]-bp[2] if bp[1]=='+' else bp[2]-bp[5]) < 1e6] # max 1mb gap -> long intron
                #todo: also check that the aligned parts have not big gap or overlap
        
        if merge:
            for i in merge:
                if new_chim[1][i][1]=='+':      #merge into next part      
                    new_chim[1][i+1][2]=new_chim[1][i][2]+new_chim[1][i+1][2]
                    new_chim[1][i+1][3]=new_chim[1][i][3]+new_chim[1][i+1][3]        
                else:
                    new_chim[1][i+1][2]=new_chim[1][i+1][2]+new_chim[1][i][2]
                    new_chim[1][i+1][3]=new_chim[1][i+1][3]+new_chim[1][i][3]
            new_chim[1]=[part for i,part in enumerate(new_chim[1]) if i not in merge]           
            bpts=tuple(bp for i,bp in enumerate(bpts)  if i not in merge)
        if bpts:
            chimeric_dict.setdefault(bpts,[]).append(new_chim)
        else:
            assert len(new_chim[1])==1
            non_chimeric.append(new_chim)
    return chimeric_dict, non_chimeric


def _add_sample_transcript(self,tr,chrom,sample_name,fuzzy_junction=5):
    'add transcript to gene in chrom - return gene on success and None if no Gene was found'
    if chrom not in self.data:
        tr['annotation']  =(4,{'intergenic':[]})
        return None
    genes_ol =self.data[chrom][tr['exons'][0][0]: tr['exons'][-1][1]]
    genes_ol_strand= [g for g in genes_ol if g.strand==tr['strand'] ] 
    #check if transcript is already there (e.g. from other sample):
    for g in genes_ol_strand:
        for tr2 in g.transcripts:
            if splice_identical(tr2['exons'], tr['exons']):
                tr2['coverage'][sample_name]=tr['coverage'][sample_name]
                if 'long_intron_chimeric' in tr:
                    tr2.setdefault('long_intron_chimeric',[]).append(tr['long_intron_chimeric'])
                #todo:check other infos that might get lost here
                return g
    #check if gene is already there (e.g. from same or other sample):
    #g,_=_get_intersects(genes_ol, tr['exons'])
    g,ref_ol,additional,not_covered=_find_splice_sites(genes_ol_strand, tr['exons'])
    if g is not None:
        if g.is_annotated: #check for fuzzy junctions (e.g. same small shift at 5' and 3' compared to reference annotation)
            shifts=g.correct_fuzzy_junctions(tr, fuzzy_junction, modify=True) #this modifies tr['exons']
            if shifts: 
                tr.setdefault('fuzzy_junction', {}).setdefault(sample_name,[]).append(shifts) #keep the info, mainly for testing/statistics    
                for tr2 in g.transcripts:#check if correction made it identical to existing
                    if splice_identical(tr2['exons'], tr['exons']):
                        tr2.setdefault('fuzzy_junction', {}).setdefault(sample_name,[]).append(shifts) #keep the info, mainly for testing/statistics    
                        tr2['coverage'][sample_name]=tr2['coverage'].get(sample_name,0)+tr['coverage'][sample_name]
                        return g 
            tr['annotation']=g.ref_splice_graph.get_alternative_splicing(tr['exons'], g.strand, additional)
            if not_covered:
                tr['novel_splice_sites']=not_covered #todo: might be changed by fuzzy junction
            #{'sj_i': sj_i, 'base_i':base_i,'category':SPLICE_CATEGORY[altsplice[1]],'subcategory':altsplice[1]} #intersects might have changed due to fuzzy junctions

        else: #add to existing novel (e.g. not in reference) gene
            start, end=min(tr['exons'][0][0],g.start),max(tr['exons'][-1][1],g.end)
            tr['annotation']=(4,_get_novel_type(genes_ol, genes_ol_strand, ref_ol))
            if start<g.start or end>g.end:# range of the novel gene might have changed 
                new_gene=Gene(start, end,g.data,self)
                self.data[chrom].add(new_gene)#todo: potential issue: in this case two genes may have grown together
                self.data[chrom].remove(g)
                g=new_gene
        #if additional:
        #    tr['annotation']=(4,tr['annotation'][1]) #fusion transcripts... todo: overrule tr['annotation']
        g.data.setdefault('transcripts',[]).append(tr) 
        g.data['splice_graph']=None #gets recomputed on next request      
        g.data['coverage'] =None 
    else:
        tr['annotation']  =(4,_get_novel_type(genes_ol, genes_ol_strand, ref_ol))
    return g

def _get_novel_type(genes_ol, genes_ol_strand, ref_ol):
    if len(ref_ol):
        return {'genic genomic':list(ref_ol.keys())}
    elif len(genes_ol_strand):
        return {'intronic':[g.name for g in genes_ol_strand]}
    elif len(genes_ol):
        return {'antisense':[g.name for g in genes_ol]}
    else:
        return {'intergenic':[]} 

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

def _find_splice_sites(genes_ol, exons):
    '''check the splice site intersect of all overlapping genes and return 
        1) the best gene, 2) exonic reference gene overlap 3) names of genes that cover additional splice sites, and 4) splice sites that are not covered. 
        For mono-exon genes return the gene with largest exonic overlap'''
    if genes_ol:
        ref_ol={g.name:g.ref_splice_graph.get_overlap(exons) for  g in genes_ol if g.is_annotated}
        if len(exons)==1: #no splice sites, but return gene with largest overlap
            ol=np.array([ref_ol[g.name] if g.is_annotated else g.splice_graph.get_overlap(exons) for  g in genes_ol ])
            best_idx=ol.argmax()
            if ol[best_idx]>0:
                return genes_ol[best_idx],ref_ol,None,None
        else:
            splice_sites=np.array([g.ref_splice_graph.find_splice_sites(exons) if g.is_annotated else g.splice_graph.find_splice_sites(exons) for  g in genes_ol ])
            sum_ol=splice_sites.sum(1)
            best_idx=sum_ol.argmax() #index of gene that covers the most splice sites
            if sum_ol[best_idx]>0:
                not_in_best=np.where(~splice_sites[best_idx])[0]
                additional=splice_sites[:,not_in_best]# additional= sites not covered by top gene
                elsefound=[(g.name,not_in_best[a]) for g,a in zip(genes_ol,additional) if a.sum()>0 ] # genes that cover additional splice sites
                notfound=(splice_sites.sum(0)==0).nonzero() #not covered splice sites
                return genes_ol[best_idx],ref_ol, elsefound, notfound
    else:
        ref_ol={}
    return None,ref_ol,None,np.zeros((len(exons)-1)*2, dtype=bool) 


def _get_intersects(genes_ol, exons):
    'calculate the intersects of all overlapping genes and return the best'  
    logger.warning('_get_intersects is depreciated')
    ## todo: in order to detect readthrough fusion genes, report all overlapping reference genes? or exon/junction wise?
    # prefer annotated genes
    intersect=[(i,g.ref_splice_graph.get_intersects(exons)) for i,g in enumerate(genes_ol) if g.is_annotated]
    if intersect:
        best_idx ,intersects= max(intersect, key=lambda x:x[1])
        if intersects[0] > 0 or (len(exons)==1 and intersects[1] > 0) :    # at least one splice overlap for multiexon genes
            return genes_ol[best_idx],intersects
    #no annotated overlap, look for novel genes
    intersect=[(i,g.splice_graph.get_intersects(exons)) for i,g in enumerate(genes_ol) if not g.is_annotated]
    if intersect:
        best_idx ,intersects= max(intersect, key=lambda x:x[1])
        if intersects[0] > 0 or (len(exons)==1 and intersects[1] > 0) :    # at least one splice overlap for multiexon genes
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
    return genes


def collapse_immune_genes(self, maxgap=300000):
    assert not self.samples, 'immune gene collapsing has to be called before adding long read samples'
    num={'IG':0, 'TR':0}
    for chrom in self.data:
        for strand in ('+','-'):
            immune={'IG':[], 'TR':[]}            
            for g in self.data[chrom]:
                if g.strand!=strand:
                    continue
                gtype=g.data['reference']['gene_type']
                if gtype[:2] in immune and gtype[-5:]=='_gene':
                    immune[gtype[:2]].append(g)
            for itype in immune:
                immune[itype].sort(key=lambda x: (x.start, x.end))            
                offset=0
                for i,g in enumerate(immune[itype]):
                    self.data[chrom].remove(g)
                    if i+1==len(immune[itype]) or g.end - immune[itype][i+1].start>maxgap:
                        ref_info={'gene_type':f'{itype}_gene', 'transcripts':[t for g in immune[itype][offset:i+1] for t in g.ref_transcripts]}
                        info={'ID':f'{itype}_locus_{num[itype]}', 'strand':strand,'chr':chrom, 'reference':ref_info}
                        start=immune[itype][offset].start
                        end=immune[itype][i].end
                        new_gene=Gene(start,end,info, self)
                        self.data[chrom].add(new_gene)
                        num[itype]+=1
                        offset=i+1
    logger.info(f'collapsed {num["IG"]} immunoglobulin loci and  {num["TR"]} T-cell receptor loci')
                        

        




def import_gff_transcripts(fn, transcriptome, chromosomes=None, gene_categories=['gene']):
    '''import transcripts from gff file (e.g. for a reference)
    returns a dict interval trees for the genes'''
    #file_size=os.path.getsize(fn) # does not help for eat 
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
            info = dict([pair.split('=', 1) for pair in ls[8].rstrip(';').split(";")]) # some gff lines end with ';' in gencode 36 
        except ValueError:
            logger.warning("GFF format error in infos (should be ; seperated key=value pairs). Skipping line:\n"+line)
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
    logger.debug('sorting exon positions...')
    for tid in exons:
        exons[tid].sort()
    missed_genes=[gid for gid in transcripts.keys() if gid not in gene_set]
    if missed_genes:        
        #logger.debug('/n'.join(gid+str(tr) for gid, tr in missed_genes.items()))
        notfound=len(missed_genes)
        found=sum((len(t) for t in genes.values()) )
        logger.warning('Missing genes! Found gene information in categories {} for {}/{} genes'.format(gene_categories,found, found+notfound))
    logger.debug('building gene data structure...')
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
    return genes

## io utility functions
def get_mutations_from_bam(bam_file,genome_file,region, min_cov=.05):
    '''not very efficient function to fetch mutations within a region from a bam file'''
    mut=dict()
    exons=[]
    n=0
    with AlignmentFile(bam_file, "rb") as align:        
        for read in align.fetch(*region):
            n+=1
            
            exons.append(junctions_from_cigar(read.cigartuples,read.reference_start))
            mutations=get_mutations(read.cigartuples, read.query_sequence,read.reference_start,read.query_qualities)
            for pos,ref,alt,qual in mutations:
                mut.setdefault(pos,{}).setdefault(alt,[0,ref,[]])
                mut[pos][alt][0]+=1
                if qual:
                    mut[pos][alt][2].append(qual)
    if min_cov<1:
        min_cov=n*min_cov

    mut={pos:v for pos,v in mut.items() if sum(v[alt][0] for alt in v)>min_cov }
    with FastaFile(genome_file) as genome_fh:
        for pos,v in mut.items():
            for alt in v.values():
                alt[1]='' if alt[1]<0 else genome_fh.fetch(read.reference_name, pos,pos+alt[1])
            for tr in exons:
                for e in tr:
                    if pos>=e[0] and pos<=e[1]:
                        mut[pos]['cov']=mut[pos].get('cov',0)+1
    return mut



def get_mutations(cigartuples, seq, ref_start,qual):
    'look up the bases affected by mutations as reported in the cigar string'
    #cigar numbers:
    #012345678
    #MIDNSHP=X
    mutations=[]
    ref_pos=ref_start
    seq_pos=0
    for cigar in cigartuples:
        if cigar[0] in (1,2,8):#I(ins), D(del) or X (missmatch): 
            ref=-cigar[1] if cigar[0] == 1 else cigar[1]
            alt_base='' if cigar[0] == 2 else seq[seq_pos:(seq_pos+cigar[1])]
            mutations.append((ref_pos,ref, alt_base,qual[seq_pos] if qual else None))
        if cigar[0] in (0, 2, 3, 7, 8):  # MDN=X -> move forward on reference
            ref_pos += cigar[1]
        if cigar[0] in (0, 1, 4,7, 8):  # MIS=X -> move forward on seq
            seq_pos += cigar[1]
    return mutations      

def aligned_part(cigartuples, is_reverse):
    "returns the interval of the trasncript that is aligned (e.g. not clipped) according to cigar. Positions are according to transcript strand"
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

def transcript_table(self, region=None, extra_columns=None,  include=None, remove=None, min_coverage=None, max_coverage=None): 
    'create a transcript table'
    if extra_columns is None:
        extra_columns=[]
    if not isinstance( extra_columns, list):
        raise ValueError('extra_columns should be provided as list')
    extracolnames={'annotation':('novelty_class', 'novelty_subclasses'),'coverage':(f'coverage_{sn}' for sn in self.samples)}

    colnames=['chr', 'transcript_start', 'transcript_end','strand', 'gene_id','gene_name' ,'transcript_nr']  + [n for w in extra_columns for n in (extracolnames[w] if w in extracolnames else (w,))]
    rows=[]
    for g,trid, tr in self.iter_transcripts(region, include,remove, min_coverage, max_coverage):                
        rows.append([g.chrom, tr['exons'][0][0], tr['exons'][-1][1], g.strand,g.id, g.name, trid]+g.get_infos(trid,extra_columns))
    df = pd.DataFrame(rows, columns=colnames)
    return(df)

def chimeric_table(self, region=None,  include=None, remove=None, star_chimeric=None, illu_len=200):
    'find chimeric genes and compile a table with relevant infos (breakpoints, coverage, ...)'
    #todo: correct handeling of three part fusion events not yet implemented
    #todo: ambiguous alignment handling not yet implemented
    
    if star_chimeric is None:
        star_chimeric=dict()
    assert isinstance(star_chimeric, dict)
    
    chim_tab=list()
    for bp,chimeric in self.chimeric.items():
        cov=tuple(sum(c.get(sa,0) for c,_ in chimeric) for sa in self.samples)
        genes=[info[4] if info[4] is not None else 'intergenic' for info in chimeric[0][1]]
        for i,bp_i in enumerate(bp):
            chim_tab.append(('_'.join(genes),)+bp_i[:3]+(genes[i],)+bp_i[3:]+(genes[i+1],)+(sum(cov),)+cov)
    chim_tab=pd.DataFrame(chim_tab,columns=['name','chr1','strand1','breakpoint1', 'gene1','chr2','strand2','breakpoint2','gene2','total_cov']+[s+'_cov' for s in self.infos['sample_table'].name])
    return chim_tab

    #todo: integrate short read coverage from star files    
    breakpoints={} #todo: this should be the isoseq breakpoints
    offset=10+len(self.infos['sample_table'])
    for sa_idx,sa in enumerate(star_chimeric):
        star_tab=pd.read_csv(star_chimeric[sa], sep='\t')
        for _, row in star_tab.iterrows():
            if row['chr_donorA'] in breakpoints and row['chr_acceptorB'] in breakpoints:
                idx1={bp.data for bp in breakpoints[row['chr_donorA']][row['brkpt_donorA']]}
                if idx1:
                    idx2={bp.data for bp in breakpoints[row['chr_acceptorB']][row['brkpt_acceptorB']]}
                    idx_ol={idx for idx,snd in idx2 if (idx, not snd) in idx1}
                    for idx in idx_ol:
                        chim_tab[idx][offset+sa_idx]+=1
                    
    chim_tab=pd.DataFrame(chim_tab, columns=['trid','len', 'gene1','part1','breakpoint1', 'gene2','part2','breakpoint2','total_cov']+[s+'_cov' for s in self.infos['sample_table'].name]+[s+"_shortread_cov" for s in star_chimeric])
    return chim_tab        

def write_gtf(self, fn, source='isotools',use_gene_name=False,  include=None, remove=None):     
    'writes the transcripts in gtf format to a file'
    
    with open(fn, 'w') as f:     
        logger.info(f'writing gtf file to {fn}')
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