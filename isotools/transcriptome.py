#!/usr/bin/python3
import os
import os.path

import itertools
import warnings
from itertools import combinations, product,tee
#From itertools receipes:
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable) #tee splits one iterator into n
    next(b, None)
    return zip(a, b)

from collections import ChainMap
import argparse
from pysam import TabixFile, AlignmentFile, FastaFile
from Bio import pairwise2
from Bio.Seq import Seq
import pandas as pd
import numpy as np
#from Bio import SeqIO, Seq, SeqRecord
import matplotlib.pyplot as plt
from intervaltree import IntervalTree, Interval
from tqdm import tqdm
import isotools.splice_graph 
from isotools.pacbio_bam import PacbioBam
import logging
import copy
import pickle
import re
#import operator as o

log=logging.getLogger(__name__)
log.setLevel(logging.INFO)
log_format=logging.Formatter('%(levelname)s: [%(asctime)s] %(name)s: %(message)s')
#log_file=logging.FileHandler('logfile.txt')
log_stream=logging.StreamHandler()
log_stream.setFormatter(log_format)
log.handlers=[]
log.addHandler(log_stream)
#a_th=.5, rtts_maxcov=10, rtts_ratio=5, min_alignment_cov=.1
default_gene_filter={'NOVEL_GENE':'all(tr["annotation"] for tr in transcripts.values())'}
default_transcript_filter={
        'CLIPPED_ALIGNMENT':'(aligned[1]-aligned[0])/source_len<.8 and "fusion" not in locals()',#more then 20% clipped
        'A_CONTENT':'downstream_A_content>.5', #more than 50% a
        'RTTS':'template_switching and any(ts[2]<=10 and ts[2]/(ts[2]+ts[3])<.2 for ts in template_switching)', #less than 10 reads and less then 20% of total reads for at least one junction
        'NONCANONICAL_SPLICING':'noncanonical_splicing',
        'NOVEL_TRANSCRIPT':'annotation is None or "splice_identical" not in annotation["as"]',
        'TRUNCATION':'truncated',
        'REFERENCE':'annotation',# and "splice_identical" in annotation',
        'UNSPLICED':'len(exons)==1','MULTIEXON':'len(exons)>1'}

CIGAR_DICT={c:i for i,c in enumerate('MIDNSHP=X')} #this is how cigar is encoded in bam

class Transcriptome:
    '''Container class for genes
    '' contains a dict of interval trees for the genes, each containing the splice graph'''
    def __new__(cls, pickle_file=None,**kwargs):
        if pickle_file is not None:
            obj=cls.load(pickle_file)
        else:
            obj=super().__new__(cls,**kwargs)
        return obj

    def __init__(self,pickle_file=None,**kwargs ):     
        if 'data' in kwargs:
            self.data,self.infos=kwargs['data'],kwargs.get('infos',dict())
            assert 'reference_file' in self.infos 
            self.make_index()
        

    @classmethod
    def from_reference(cls, reference_file, file_format='auto',**kwargs):
        tr = cls.__new__(cls)         
        if file_format=='auto':        
            file_format=os.path.splitext(reference_file)[1].lstrip('.')
            if file_format=='gz':
                file_format=os.path.splitext(reference_file[:-3])[1].lstrip('.')
        log.info(f'importing reference from {file_format} file {reference_file}')
        if file_format == 'gtf':
            tr.data,tr.infos= import_gtf_transcripts(reference_file,tr,  **kwargs)
        elif file_format in ('gff', 'gff3'):
            tr.data,tr.infos= import_gff_transcripts(reference_file,tr,  **kwargs)
        elif file_format == 'pkl':
            genes,infos= pickle.load(open(reference_file, 'rb'))
            tr=next(iter(next(iter(genes.values()))))._transcriptome            
            if [k for k in infos if k!='reference_file']:
                log.warning('the pickle file seems to contain additional expression information... extracting refrence')
                ref_data=tr._extract_reference()
                tr.data,tr.infos=tr._extract_reference(), {'reference_file':infos['reference_file']}
        return tr

    @classmethod
    def load(cls, pickle_file):
        'restores the information of a transcriptome from a pickle file'
        log.info('loading transcriptome from '+pickle_file)
        data, infos=pickle.load(open(pickle_file, 'rb'))
        tr=next(iter(next(iter(data.values()))))._transcriptome
        return tr

        
    def save_reference(self, fn=None):    
        'saves the reference information of a transcriptome in a pickle file'
        if fn is None:
            fn=self.infos['reference_file']+'.isotools.pkl'
        log.info('saving reference to '+fn)       
        ref_data=self._extract_reference() 
        pickle.dump((ref_data,{'reference_file':self.infos['reference_file']} ))

    def _extract_reference(self):
        if not [k for k in self.infos if k!='reference_file']:
            return self.data #do nothing
        ref_data={} # extract the reference
        for chrom,tree in self.data.items():
            ref_data[chrom]=IntervalTree(Gene(g.start,g.end,{k:g.data[k] for k in Gene.required_infos+['reference']}, self) for g in tree if g.is_annotated)
        return ref_data

    def save(self, fn=None):
        'saves the information of a transcriptome (including reference) in a pickle file'
        if fn is None:
            fn=self.infos['file_name']+'.isotools.pkl'
        log.info('saving transcriptome to '+fn)
        pickle.dump((self.data,self.infos), open(fn, 'wb'))
    
    
    def write_gtf(self, fn, source='isotools',use_gene_name=False,  include=None, remove=None):     
        'writes the transcripts in gtf format to a file'
        with open(fn, 'w') as f:     
            for gene in tqdm(self):
                lines=gene.to_gtf(source=source,use_gene_name=use_gene_name, include=include, remove=remove)
                if lines:
                    _=f.write('\n'.join( ('\t'.join(str(field) for field in line) for line in lines) )+'\n')

    def make_index(self):
        'updates the index used for __getitem__, e.g. the [] operator'
        idx=dict()
        for g in self:
            if g.id in idx: # at least id should be unique - maybe raise exception?
                log.warn(f'{g.id} seems to be ambigous: {str(self[g.id])} vs {str(g)}')
            idx[g.name] = g
            idx[g.id]=g
        self._idx=idx
        
    def __getitem__(self, key):
        return self._idx[key]

    def __len__(self):
        return self.n_genes
    
    def __contains__(self, key):
        return key in self._idx
    
    def remove_chromosome(self, chromosome):
        'deletes the chromosome from the transcriptome'
        del self.data[chromosome]
        self.make_index()
    
    def add_illumina_coverage(self, illumina_fn, load=False):
        'Add illumina coverage to the genes.\n This does, by default, not actually read the bams, but reading is done at first access'
        log.warning('illumina coverage untested')
        if load: # when loading coverage for all genes keep the filehandle open, hopefully a bit faster
            for bamfile in illumina_fn.values():            
                log.info(bamfile)
                with AlignmentFile(bamfile, "rb") as align:
                    try:
                        for g in tqdm(self):
                            g.data.setdefault('illumina',list()).append(Coverage.from_alignment(align,g))
                    except:
                        for g in self:
                            g.data.pop('illumina', None)
                        raise
        else:
            for g in self:
                g.data['illumina']=[Coverage.from_bam(bamfile,g, load=load) for bamfile in illumina_fn.values()]
        self.infos['illumina_fn']=illumina_fn


    def add_biases(self, genome_fn):
        'populates transcript["biases"] information, which can be used do create filters'
        with FastaFile(genome_fn) as genome_fh:
            for g in tqdm(self):                
                ts_candidates=g.splice_graph.find_ts_candidates()
                for start, end, js, ls, idx in ts_candidates:
                    for tr in (g.transcripts[i] for i in idx):
                        tr.setdefault('template_switching',[]).append((start, end, js, ls))
                g.add_direct_repeat_len(genome_fh) 
                g.add_noncanonical_splicing(genome_fh)
                g.add_threeprime_a_content(genome_fh)
                g.add_truncations()
        self.infos['biases']=True # flag to check that the function was called

    def add_filter(self, transcript_filter={},gene_filter={}):
        'create filter flags which can be used by iter_transcripts'
        #possible filter flags: 'A_CONTENT','RTTS','NONCANONICAL_SPLICING','NOVEL_GENE','NOVEL_TRANSCRIPT','TRUNCATION'
        # a_th=.5, rtts_maxcov=10, rtts_ratio=5, min_alignment_cov=.1
        #biases are evaluated with respect to specific thresholds
        gene_attributes={k for g in self for k in g.data.keys() }
        tr_attributes={k for g in self for tr in g.transcripts for k in tr.keys() }
        tr_attributes.add('filter')
        if not gene_filter and not transcript_filter:
            gene_filter=default_gene_filter
            transcript_filter=default_transcript_filter
        gene_ffun={label:filter_function(gene_attributes, fun) for label,fun in gene_filter.items()}
        tr_ffun={label:filter_function(tr_attributes, fun) for label,fun in transcript_filter.items()}
        for g in tqdm(self):
            g.add_filter(tr_ffun,gene_ffun)
        self.infos['filter']={'gene_filter':gene_filter, 'transcript_filter':transcript_filter}


    def get_sample_idx(self, group_column='name'):
        'returns a dict with group names as keys and index lists as values'
        return self.infos['sample_table'].groupby(group_column).groups

    @property
    def sample_table(self):
        try:
           return self.infos['sample_table']
        except KeyError:
            return pd.DataFrame(columns=['name','file','group'])
    
    @property
    def samples(self):
        return list(self.sample_table.name)

    @property
    def groups(self):
        return dict(self.sample_table.groupby('group')['name'].apply(list))

        
    @property
    def n_transcripts(self):
        if self.data==None:
            return 0
        return sum(g.n_transcripts for g in self)

    

    @property
    def n_genes(self):
        if self.data==None:
            return 0
        return sum((len(t) for t in self.data.values()))
    
    @property
    def novel_genes(self): #this is used for id assignment
        try:
            return self.infos['novel_counter']
        except KeyError:
            self.infos['novel_counter']=0
            return 0

    @property
    def chromosomes(self):
        return list(self.data)            

    def __str__(self):
        return '{} object with {} genes and {} transcripts'.format(type(self).__name__, self.n_genes, self.n_transcripts)
    
    def __repr__(self):
        return object.__repr__(self)

    def __iter__(self):
        return (gene for tree in self.data.values() for gene in tree)

       
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

    def add_sample_from_bam(self,fn, sample_name,fuzzy_junction=5, **kwargs):
        '''import expressed transcripts from bam and add it to existing transcriptome'''
        #todo: one alignment may contain several samples - this is not supported at the moment
        assert sample_name not in self.samples, 'sample %s is already in the data set.' % sample_name
        log.info(f'adding sample {sample_name}')
        kwargs['name']=sample_name
        kwargs['file']=fn
        self.infos['sample_table']=self.sample_table.append(kwargs, ignore_index=True)
        chromosomes=self.chromosomes
        with AlignmentFile(fn, "rb") as align:        
            stats = align.get_index_statistics()
            # try catch if sam/ no index /not pacbio?
            if not chromosomes:
                chromosomes=align.references
            total_reads = sum([s.mapped for s in stats if s.contig in chromosomes])
            
            #runs=set()
            #genes={c:IntervalTree() for c in chromosomes}
            chimeric=dict()
            
            with tqdm(total=total_reads, unit='transcripts') as pbar:
                for chrom in chromosomes: #what about secondary alignments to non listed chromosomes?
                    pbar.set_postfix(chr=chrom)
                    transcripts=IntervalTree()
                    novel=IntervalTree()
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
                        if 'SA' in tags:#part of a chimeric alignment
                            chimeric.setdefault(read.query_name, []).append({'chr':chrom,'strand':strand,'range':tr_range,'cigar':read.cigarstring, 'aligned': aligned_part(read.cigartuples, read.is_reverse)}) #cigar is required to sort out the order 
                            tr['samples'][sample_name].setdefault('chimeric',[]).append(chimeric[read.query_name])                       
                        elif 4 in read.cigartuples: #clipping
                            clip=get_clipping(read.cigartuples, read, read.reference_start)
                            tr['samples'][sample_name].setdefault('clipping',{}).setdefault(clip,0)
                            tr['samples'][sample_name]['clipping'][clip]+=cov
                    for tr_interval in transcripts:
                        tr=tr_interval.data
                        tr['exons'][0][0]=min(r[0] for r in tr['samples'][sample_name]['range'])
                        tr['exons'][-1][1]=max(r[1] for r in tr['samples'][sample_name]['range'])  
                        gene=self._add_sample_transcript(tr,chrom, sample_name) #tr is not updated
                        if gene is None:
                            novel.add(tr_interval)    
                    self._add_novel_genes(novel,chrom)
        for g in self:
            if 'splice_graph' in g.data and g.data['splice_graph'] is not None: # still valid splice graphs no new transcripts - add a row of zeros to weights
                n_samples=len(self.samples)
                w=g.coverage
                if w.shape[0]<n_samples:                    
                    w_new=np.zeros((n_samples,w.shape[1]))
                    w_new[:w.shape[0],:] = w
                    g.data['splice_graph'].weights=w_new
    
    def _add_sample_transcript(self,tr,chrom,sample_name,fuzzy_junction=5):
        'add transcript to gene in chrom - return True on success and False if no Gene was found'
        genes_ol = [g for g in self.data[chrom][tr['exons'][0][0]: tr['exons'][-1][1]] if g.strand==tr['strand'] ]
        #check if transcript is already there (e.g. from other sample):
        for g in genes_ol:
            for tr2 in g.transcripts:
                if splice_identical(tr2['exons'], tr['exons']):
                    tr2['samples'][sample_name]=tr['samples'][sample_name]
                    return g
        #check if gene is already there (e.g. from same or other sample):
        g,(sj_i, base_i)=self.get_intersects(genes_ol, tr['exons'], tr['strand'])
        if g is not None:            
            if g.is_annotated:                
                corrected=g.correct_fuzzy_junctions(tr,g, fuzzy_junction)
                if corrected: #check if correction made it identical to existing
                    return g
                altsplice= g.ref_splice_graph.get_alternative_splicing(tr['exons'], g.strand)
                tr['annotation']={'sj_i': sj_i, 'base_i':base_i,'as':altsplice}

            else: # range of the novel gene might have changed
                start, end=min(tr['exons'][0][0],g.start),max(tr['exons'][-1][1],g.end)
                if start<g.start or end>g.end:
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
            self.data[chrom].add(Gene(start,end,new_data,self ))
            self.infos['novel_counter']=n_novel
                            
    def get_intersects(self,genes_ol, exons, strand):
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
                    
'''
    def annotate(self, reference,fuzzy_junction=5,novel_params=None): 
        'depreciated, broken, replaced by new add_sample function... kept as not all functionality transfered yet'
        'Add annotation to the transcripts and alines the junctions to correct for ambigious alignments' 
        raise ValueError('depreciated')
        genes=dict()
        unknown_fusions=list()
        novel_default=dict(spj_iou_th=0, reg_iou_th=0, gene_prefix='PB_novel_')
        if novel_params is not None:
            novel_default.update(novel_params)
        offset=0
        with tqdm(total=len(self), unit='transcripts') as pbar:
            for chrom, tree in self.data.items():
                gene_idx=dict()
                novel=IntervalTree()
                pbar.set_postfix(chr=chrom)
                for g in tree:    
                    for trid,tr in g.transcripts.items():    
                        log.debug('annotateing transcript %s',trid)
                        ref_g,(sj_i, base_i)=get_annotation(chrom, tr['exons'],g.strand, reference)
                        if ref_g is not None:
                            g.correct_fuzzy_junctions(trid,ref_g, fuzzy_junction)
                            altsplice= ref_g.splice_graph.get_alternative_splicing(tr['exons'], g.strand)
                            tr['annotation']={'ref_gene_name': ref_g.name, 'ref_gene_id':ref_g.id, 'sj_i': sj_i, 'base_i':base_i,'as':altsplice}
                            gene_idx.setdefault(ref_g.name,{'chr':chrom, 'ID':ref_g.id, 'Name': ref_g.name, 'strand':g.strand ,'transcripts': dict()})
                            gene_idx[ref_g.name]['transcripts'][trid]=tr
                        else:
                            tr['annotation']=None
                        if 'fusion' in tr:
                            log.debug('annotateing fusion of %s',trid)
                            for f in tr['fusion']:
                                f_ref_g,(f_sj_i, f_base_i)=get_annotation(f['chr'],f['exons'],f['strand'],reference )
                                if f_ref_g is not None:
                                    f['annotation']={'ref_gene_name':f_ref_g.name,'ref_gene_id':f_ref_g.id, 'sj_i':f_sj_i,'base_i':f_base_i}
                                else:
                                    unknown_fusions.append(f)
                    if all(tr['annotation'] is None for tr in g.transcripts.values()):
                        novel.add(g)
                    pbar.update(1)
                    
                #merge the novel genes for the chromosome
                genes[chrom]=collapse_transcripts_to_genes(novel,offset=offset,**novel_default)
                offset+=len(genes[chrom])
                #add the genes with annotation
                for gname, g in gene_idx.items():
                    start=min(t['exons'][0][0] for t in g['transcripts'].values())
                    end=max(t['exons'][-1][1] for t in g['transcripts'].values())
                    genes[chrom].add(Gene(start, end, g,self))        
        self.data=genes
        for f in unknown_fusions:
            f_ref_g,(f_sj_i, f_base_i)=get_annotation(f['chr'],f['exons'],f['strand'],self )
            if f_ref_g is not None:
                f['annotation']={'ref_gene_name':f_ref_g.name,'ref_gene_id':f_ref_g.id, 'sj_i':f_sj_i,'base_i':f_base_i}
            else:
                f['annotation']=None
        self.make_index()
        self.infos['annotation']=reference.infos['file_name']
'''   

    def iter_genes(self, region=None):
        'iterate over the genes of a region'
        if region is None:
            genes=self
        elif region in self.data:
            genes=self.data[region] #provide chromosome
        else:
            try:
                chrom, pos=region.split(':')
                start, end=[int(v) for v in pos.split('-')]
            except ValueError:
                chrom,start,end=region
            except:
                raise ValueError('incorrect region {} - specify as string "chr:start-end" or tuple ("chr",start,end)'.format(region))
            else:
                genes=self.data[chrom][start:end]
        for g in genes:
            #todo: add filtering?
            yield g
    
    def iter_transcripts(self,region=None,include=None, remove=None):
        'iterate over the transcripts of a region, optionally applying filters'
        for g in self.iter_genes(region):
            for i,tr in g.filter_transcripts(include, remove):
                yield g,i,tr


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
        log.warn('fusion table is broken')
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

class Coverage: 
    'stores the illumina read coverage of a gene'
    #plan: make a binned version, or use run length encoding
    def __init__(self, cov, junctions, offset, chrom=None):
        self._cov=cov
        self._junctions=junctions
        self.reg=None if cov is None else (chrom,offset, offset+len(cov)) 
        self.bam_fn=None
    
    @classmethod
    def from_bam(cls, bam_fn,g, load=False):
        'assign the bam file'
        if load:
            with AlignmentFile(bam_fn, 'rb') as align:
                return cls.from_alignment(g,align)
        else: #load on demand
            obj = cls.__new__(cls)  
            #obj.__init__(None, None, None)
            obj._cov=None
            obj._junctions=None
            obj.bam_fn=bam_fn
            obj.reg=(g.chrom, g.start,g.end)
            return obj

    @classmethod
    def from_alignment(cls, align_fh,g):
        'load the coverage from bam file'
        cov,junctions=cls._import_coverage(align_fh,(g.chrom, g.start,g.end))
        obj = cls.__new__(cls)  
        obj.__init__(cov,junctions, g.start )
        return obj

    @classmethod #this is slow - called only if coverage is requested
    def _import_coverage(cls,align_fh, reg):
        delta=np.zeros(reg[2]-reg[1])        
        junctions={}
        for read in align_fh.fetch(*reg):        
            exons=junctions_from_cigar(read.cigartuples,read.reference_start) 
            #alternative: read.get_blocks() should be more efficient.. -todo: is it different?
            for i, exon in enumerate(exons):
                s=max(reg[1],min(reg[2]-1,exon[0]))-reg[1]
                e=max(reg[1],min(reg[2]-1,exon[1]))-reg[1]
                delta[s]+=1
                delta[e]-=1            
                if i>0:
                    jpos=(exons[i-1][1],exon[0])
                    if jpos[1]-jpos[0]<1:
                        continue
                    junctions[jpos]=junctions.get(jpos,0)+1            
        cov=np.cumsum(delta) #todo: use rle instead?  
        return cov,junctions
    
    def load(self):
        'load the coverage from bam file'
        with AlignmentFile(self.bam_fn, 'rb') as align:
            log.debug(f'Illumina coverage of region {self.reg[0]}:{self.reg[1]}-{self.reg[2]} is loaded from {self.bam_fn}') #info or debug?
            self._cov, self._junctions=type(self)._import_coverage(align, self.reg)
            
    @property
    def junctions(self):
        if self._junctions is None:
            self.load()
        return self._junctions
    
    @property
    def profile(self):
        if self._cov is None:
            self.load()
        return self._cov #todo: implement rle?

    def __getitem__(self, subscript): 
        if isinstance(subscript, slice): 
            return self.profile[slice(None if subscript.start is None else subscript.start-self.reg[1],
                                        None if subscript.stop is None else subscript.stop-self.reg[1],
                                        subscript.step)] #does not get extended if outside range       
        elif subscript < self.reg[1] or subscript >= self.reg[2]:
            log.warning('requested coverage outside range')
            return None
        else:
            return(self.profile[subscript-self.reg[1]])
        
            
        

class Gene(Interval):
    'this class stores all gene information and transcripts. It is derived from intervaltree.Interval'
    required_infos=['ID','name', 'chr', 'strand']

    def __new__(cls,begin, end, data, transcriptome):
        return super().__new__(cls,begin, end, data) #required as Interval (and Gene) are immutable

    def __init__(self,begin, end, data, transcriptome):
        self._transcriptome=transcriptome


    def __str__(self):
        return('Gene {} {}({}), {} reference transcripts, {} expressed transcripts'.format(self.name, self.region, self.strand, self.n_ref_transcripts, self.n_transcripts))
    
    def __repr__(self):
        return object.__repr__(self)

    def add_filter(self, transcript_filter,gene_filter):   
        'add the filter flags'
        gene_tags=[name for name,fun in gene_filter.items() if  fun(**self.data)]
        for tr in self.transcripts:
            tr['filter']=gene_tags.copy() #todo: introduce gene filter flags (in addition to the transcript filter)
            tr['filter'].extend([name for name,fun in transcript_filter.items() if fun(**tr)])

    def filter_transcripts(self, include=None, remove=None):
        ''' iterator over the transcripts of the gene. Filtering implemented by passing lists of flags to the parametrs:
        include: transcripts must have at least one of the flags
        remove: transcripts must not have one of the flags (higher priority)'''
        for i,tr in enumerate(self.transcripts):
            select =  include is None or any(f in tr['filter'] for f in include)
            select &= remove is None or all(f not in tr['filter'] for f in remove)
            if select:
                yield i,tr
        
    def correct_fuzzy_junctions(self,tr,ref_g, size, modify=True):
        'this corrects for splicing shifts (e.g. same difference compared to reference annotaiton at both donor and acceptor) presumably caused by ambigous alignments'
        exons=tr['exons']
        shifts=self.ref_splice_graph.fuzzy_junction(exons,size)
        if shifts and modify:
            for i,sh in shifts.items():
                exons[i][1]+=sh
                exons[i+1][0]+=sh
                assert len(tr['samples'])==1 
                sample_name=next(iter(tr['samples']))
                tr['samples']['fuzzy_junction']=shifts #keep the info, mainly for testing/statistics    
                for tr2 in self.transcripts:
                    if splice_identical(tr2['exons'], exons):
                        tr2['samples'][sample_name]=tr['samples'][sample_name]
        return shifts
        

    def to_gtf(self, source='isoseq', use_gene_name=False, include=None, remove=None):
        'creates the gtf lines of the gene as strings'
        info={'gene_id':self.name if use_gene_name else self.id}
        gene=(self.chrom, source, 'gene', self.start+1, self.end, '.',self.strand, '.', '; '.join('{} "{}"'.format(k,v) for k,v in info.items() ))
        lines=list()
        for i,tr in enumerate(self.filter_transcripts(include, remove)):
            info['transcript_id']=f'{self.id}.{i}'
            #todo: print only relevant infos
            lines.append((self.chrom, source, 'transcript', tr['exons'][0][0]+1, tr['exons'][-1][1], '.',self.strand, '.', '; '.join('{} "{}"'.format(k,v) for k,v in info.items() if k != 'exon_id')))
            for enr, pos in enumerate(tr['exons']):
                info['exon_id']='{}.{}_{}'.format(self.id,i, enr)
                lines.append((self.chrom, source, 'exon', pos[0]+1, pos[1], '.',self.strand, '.', '; '.join('{} "{}"'.format(k,v) for k,v in info.items())))
        if lines:
            lines.append(gene)
        return(lines)
    


    def add_noncanonical_splicing(self, genome_fh):
        '''for all transcripts, add the the index and sequence of noncanonical (i.e. not GT-AG) splice sites
        i.e. save the dinucleotides of donor and aceptor as XX-YY string in the transcript biases dicts
        noncanonical splicing might indicate technical artifacts (template switching, missalignment, ...)'''
        
        for tr in self.transcripts:
            pos=((tr['exons'][i][1], tr['exons'][i+1][0]-2) for i in range(len(tr['exons'])-1))
            sj_seq=((genome_fh.fetch(self.chrom, p, p+2).upper() for p in i) for i in pos)
            if self.strand=='+':
                sj_seq=[d+a for d,a in sj_seq]
            else:
                sj_seq=[str(Seq(d+a).reverse_complement()) for d,a in sj_seq]
            nc=[(i,seq) for i,seq in enumerate(sj_seq) if seq != 'GTAG']
            if nc:
                tr['noncanonical_splicing']=nc

    def add_direct_repeat_len(self, genome_fh,delta=10):
        'perform a global alignment of the sequences at splice sites and report the score. Direct repeats indicate template switching'
        for tr in self.transcripts:
            introns=((tr['exons'][i][1], tr['exons'][i+1][0]) for i in range(len(tr['exons'])-1))
            intron_seq=((genome_fh.fetch(self.chrom, pos-delta, pos+delta) for pos in i) for i in introns)
            #intron_seq=[[genome_fh.fetch(self.chrom, pos-delta, pos+delta) for pos in i] for i in introns]
            #align=[pairwise2.align.globalms(seq5, seq3, 1,-1, -1, 0) for seq5, seq3 in intron_seq] 
            #align=[pairwise2.align.globalxs(seq5, seq3,-1,-1, score_only=True) for seq5, seq3 in intron_seq] # number of equal bases at splice site
            align=[[a==b for a,b in zip(*seq)] for seq in intron_seq ]
            score=[0]*len(align)
            for i,a in enumerate(align): #get the length of the direct repeat at the junction (no missmatches or gaps)
                pos=delta
                while(a[pos]):
                    score[i]+=1
                    pos+=1
                    if pos==len(a):
                        break
                pos=delta-1
                while(a[pos]):
                    score[i]+=1
                    pos-=1
                    if pos<0:
                        break
            #parameters 2,-3,-1,-1 = match , missmatch, gap open, gap extend
            #this is not optimal, since it does not account for the special requirements here: 
            #a direct repeat shoud start at the same position in the seqs, so start and end missmatches free, but gaps not
            #tr['ts_score']=list(zip(score,intron_seq))
            tr['direct_repeat_len']=[min(s, delta) for s in score]
    
    def add_threeprime_a_content(self, genome_fh, length=30):
        'add the information of the genomic A content downstream the transcript. High values indicate genomic origin of the pacbio read'
        for tr in self.transcripts:
            if self.strand=='+':
                pos=tr['exons'][-1][1]
            else:
                pos=tr['exons'][0][0]-length
            seq=genome_fh.fetch(self.chrom, max(0,pos), pos+length) 
            if self.strand=='+':
                tr['downstream_A_content']=seq.upper().count('A')/length
            else:
                tr['downstream_A_content']=seq.upper().count('T')/length

    def add_truncations(self): 
        'check for truncations and save indices of truncated genes in transcripts'
        #dist5=-1, dist3=10
        is_truncated=set()
        for trid1, trid2 in combinations(range(len(self.transcripts)), 2):
            if trid1 in is_truncated or trid2 in is_truncated:
                continue
            if len(self.transcripts[trid1]['exons'])>len(self.transcripts[trid2]['exons']):
                trid1, trid2=trid2, trid1
            short_tr=self.transcripts[trid1]['exons']
            long_tr=self.transcripts[trid2]['exons']
            truncation, delta1, delta2=is_truncation(long_tr, short_tr)
            if truncation:
                #handle splice identical cases, test if trid1 is the long version 
                if delta1+delta2<0:
                    trid1, trid2=trid2, trid1
                    delta1*=-1
                    delta2*=-1
                is_truncated.add(trid1)                
                if self.strand=='+':
                    #if (dist5<0 or dist5>abs(delta1)) and (dist3<0 or dist3>abs(delta2)):
                    self.transcripts[trid1]['truncated']=(trid2,delta1, delta2)
                else:
                    #if (dist5<0 or dist5>abs(delta2)) and (dist3<0 or dist3>abs(delta1)):
                    self.transcripts[trid1]['truncated']=(trid2,delta2, delta1)
        
    def coding_len(self, trid):
        'returns length of 5\'UTR, coding sequence and 3\'UTR'
        try:            
            exons=self.transcripts[trid]['exons']
            cds=self.transcripts[trid]['CDS']
        except KeyError:
            return None
        else:
            coding_len=_coding_len(exons, cds)
        if self.strand=='-':
            coding_len.reverse()
        return coding_len

    def get_values(self ,trid, what):
        'get the gene information specified in "what" as a tuple'
        if what is None:
            return ()
        return tuple((v for w in what for v in self._get_value(trid, w)))

    def _get_value(self, trid, what):
        if what=='length':
            return sum((e-b for b,e in self.transcripts[trid]['exons'])),#return a tuple (hence the comma)
        if what=='filter':
            if self.transcripts[trid]['filter']:
                return ','.join(self.transcripts[trid]['filter']),
            else:
                return 'PASS',
        elif what=='n_exons':
            return len(self.transcripts[trid]['exons']),
        elif what=='exon_starts':
            return ','.join(str(e[0]) for e in self.transcripts[trid]['exons']),
        elif what=='exon_ends':
            return ','.join(str(e[1]) for e in self.transcripts[trid]['exons']),
        elif what=='annotation':
            #sel=['sj_i','base_i', 'as']
            support=self.transcripts[trid]['annotation']
            if support is None:
                return ('NA',)*3
            else:
                #vals=support[n] if n in support else 'NA' for n in sel
                stype=support['as']
                if isinstance(stype, dict):
                    type_string=';'.join(k if v is None else '{}:{}'.format(k,v) for k,v in stype.items())
                else:
                    type_string=';'.join(str(x) for x in stype)
                vals=(support['sj_i'],support['base_i'],type_string)
                return(vals)
        elif what=='coverage':
            return self.transcripts[trid]['coverage']
        elif what=='downstream_A_content':
            return self.transcripts[trid]['downstream_A_content'],
        elif what in self.transcripts[trid]:
            val=str(self.transcripts[trid][what])
            try:
                iter(val)
            except TypeError:
                return val, #atomic (e.g. numeric)
            else:
                return str(val), #iterables get converted to string
        return '',

    @property
    def coverage(self):
        return self.splice_graph.weights


    @property
    def gene_coverage(self):
        return self.coverage.sum(1)

    @property
    def chrom(self):
        try:
            return self.data['chr']
        except KeyError: 
            raise
 
    @property
    def start(self): #alias for begin
        return self.begin
    
    @property
    def region(self):
        try:
            return '{}:{}-{}'.format(self.chrom,self.start, self.end )
        except KeyError:
            raise

    @property 
    def illumina_coverage(self):
        try:
            return self.data['illumina']
        except KeyError:
            raise ValueError(f'no illumina coverage for gene {self.name} add illumina bam files first')

    @property
    def id(self):
        try:
            return self.data['ID']    
        except KeyError:
            log.error(self.data)
            raise
            
    @property
    def name(self):
        try:
            return self.data['name']
        except KeyError:
            return self.id #e.g. novel genes do not have a name (but id)
    
    @property
    def is_annotated(self):
        return 'reference' in self.data
    
    @property
    def is_expressed(self):
        return bool(self.transcripts)


    @property
    def strand(self):
        return self.data['strand']
        

    @property
    def transcripts(self):
        try:
            return self.data['transcripts']
        except KeyError:
            return []

    @property
    def ref_transcripts(self):
        try:
            return self.data['reference']['transcripts']
        except KeyError:
            return {}

    @property
    def n_transcripts(self):
        return len(self.transcripts)
    
    @property
    def n_ref_transcripts(self):
        return len(self.ref_transcripts)

    @property
    def ref_splice_graph(self): #raises key error if not self.is_annotated
        assert self.is_annotated, "reference splice graph requested on novel gene"
        if 'splice_graph' not in self.data['reference'] or self.data['reference']['splice_graph'] is None:
            exons=[tr['exons'] for tr in self.ref_transcripts]
            self.data['reference']['splice_graph']=isotools.splice_graph.SpliceGraph(exons)
        return self.data['reference']['splice_graph']
        
    @property
    def splice_graph(self):
        if 'splice_graph' not in self.data or self.data['splice_graph'] is None:
            self.data['splice_graph']=isotools.splice_graph.SpliceGraph(self, samples=self._transcriptome.samples)
        return self.data['splice_graph']

    def __copy__(self):
        return Gene(self.start, self.end, self.data, self._transcriptome)        
        

    def __deepcopy__(self, memo):
        return Gene(self.start, self.end, copy.deepcopy(self.data, memo), self._transcriptome)
        

    def __reduce__(self):
        return Gene, (self.start, self.end, self.data, self._transcriptome)  

    def copy(self):
        return self.__copy__()

def genomic_position(exons, tr_pos):
    'gets a sorted list of transcript positions and finds the corresponding genomic positions'
    offset=0
    g_pos=list()
    exon_iter=iter(exons)
    for pos in tr_pos:
        while offset<pos:
            try:
                e=next(exon_iter)        
            except StopIteration:
                raise ValueError(f'tr_pos {pos} outside transcript of length {sum(e[1]-e[0] for e in exons)}')
            offset+=e[1]-e[0]
        g_pos.append(e[1]+pos-offset)
    return g_pos

def _coding_len(exons, cds):
    coding_len=[0,0,0]            
    state=0
    for e in exons:
        if state<2 and e[1]>=cds[state]:
            coding_len[state]+=cds[state]-e[0]
            if state==0 and cds[1]<=e[1]: # special case: CDS start and end in same exon
                coding_len[1]=cds[1]-cds[0]
                coding_len[2]=e[1]-cds[1]
                state+=2
            else:
                coding_len[state+1]=e[1]-cds[state]
                state+=1
        else:
            coding_len[state]+=e[1]-e[0]
    return coding_len


#io
def query_bam(bam_fn,use_pbi=True,**kwargs):
        with PacbioBam.open(bam_fn, use_pbi) as align:
            for read in align.fetch(**kwargs):
                yield read
                

 

def import_pacbio_transcripts(fn,sample_table,genome_fn=None, chromosomes=None, orf=False):
    '''import transcripts from pacbio bam
    'returns a dict interval trees for the genes, each containing the splice graph'''
    raise Exception('import_pacbio_transcripts is broken/depreciated')
    n_secondary=0
    neglected_runs=set()
    #runs=pd.DataFrame(runs)
    assert isinstance(sample_table, pd.DataFrame), '\"sample_table\" must be a pandas DataFrame'
    sample_table=sample_table.reset_index(drop=sample_table.index.name is None)# make sure the index is 0..n
    assert 'name' in sample_table.columns and 'run' in sample_table.columns, '\"sample_table\" must be a pandas DataFrame with columns \"name\" and \"run\"'
    
    run_idx={v:k for k,v in sample_table.run.items()}
    #with FastaFile(genome_fn) if genome_fn is not None else None as genome_fh: #none has no open, does not work
    try:
        genome_fh=FastaFile(genome_fn) if genome_fn is not None else None
        with AlignmentFile(fn, "rb") as align:        
            stats = align.get_index_statistics()
            # try catch if sam/ no index /not pacbio?
            
            if chromosomes is None:
                chromosomes=align.references
            total_reads = sum([s.mapped for s in stats if s.contig in chromosomes])
            #runs=set()
            genes={c:IntervalTree() for c in chromosomes}
            fusion=dict()
            fusion_primary=dict()
            for read in tqdm(align, total=total_reads, unit='transcripts'):
                if read.flag & 0x100:
                    n_secondary+=1
                    continue # use only primary alignments
                chrom = read.reference_name
                if chrom not in chromosomes and not read.flag & 0x800: #allow secondary alignments form other chromosomes
                    continue
                strand = '-' if read.is_reverse else '+'
                exons = junctions_from_cigar(read.cigartuples,read.reference_start)
                tags= dict(read.tags)
                cov=np.zeros(len(sample_table))
                for read_id in tags['im'].split(','):
                    run=read_id[:read_id.find('/')]
                    try:
                        cov[ run_idx[run] ]+=1
                    except KeyError:
                        neglected_runs.add(run)              
                
                #transcript infos
                tr={'exons':exons, 'source_len':len(read.query_sequence), 'cigar':read.cigarstring, 
                    'aligned': aligned_part(read.cigartuples, read.is_reverse),
                    'coverage':cov} 
                if orf: #this takes time, hence made optional
                    tr['CDS']=find_orf(read.query_sequence)
                if genome_fh is not None:
                    tr['mutations']=get_mutations(read.cigartuples, read.query_sequence, genome_fh, chrom,read.reference_start)
                #gene infos
                info={'strand':strand, 'chr':chrom,
                        'ID':read.query_name,'transcripts':{read.query_name:tr}}
                
                if 'SA' in tags: #part of a fusion gene   
                    if read.flag & 0x800: #secondary part
                        fusion.setdefault(read.query_name, []).append({**{k:tr[k] for k in ['exons','aligned', 'cigar', 'mutations'] if k in tr},**{'chr':chrom,'strand':strand}}) #cigar is required to sort out the order 
                    else: #primary part
                        fusion_primary[read.query_name]=(tr,chrom,strand)
                    #chimeric transcripts are added to genes later, since the start and end  are not known yet
                else:
                    genes[chrom].add(Gene(exons[0][0],exons[-1][1],info, transcriptome))
    except:
        if genome_fh is not None:
            genome_fh.close()
        raise
    if genome_fh is not None:
        genome_fh.close()
    if neglected_runs:
        log.warning(f'found the follwing runs in pacbio file, which where not specified in the run table: {neglected_runs}')    
    #link secondary alignments (chimeric genes)
    for tid,(primary,chrom,strand) in fusion_primary.items():
        # distinguish fusion and long introns (introns > 200kb seem to be stored as chimeric alignments with minimap2 by default (-G))
        # todo: true fusions on same chr (e.g. readthroughs) will be detected during annotation step
        primary['exons'],primary['cigar'],primary['aligned'], kept_fusions=combine_chimeric(primary,chrom,strand,fusion[tid])
        if kept_fusions:
            primary['fusion']=kept_fusions
        info={'strand':strand, 'chr':chrom,
                        'ID':tid,'transcripts':{tid:primary}}
        genes[chrom].add(Gene(exons[0][0],exons[-1][1],info))
    if n_secondary>0:
        log.info(f"skipped {n_secondary} secondary transcripts {n_secondary/total_reads*100}%")
    return genes,{'sample_table':sample_table}

def combine_chimeric(tr,chrom,strand,fusions, tr_dist=10, max_intron=1e6):
    'resolve long introns, which are reported as chimeric alignments by the alignment tool but can be explained with splicing'
    #todo: think about muting the original tr dict (instad of return)
    exons=tr['exons'].copy()
    cigar=tr['cigar'] #no copy needed here as string is imutable
    #aligned is according to transcript strand, but we need it according to genome fwd strand to be compatible with exons
    if strand=='+':
        aligned=tr['aligned'].copy()
        f_aligned=[f['aligned'] for f in fusions] #we care about f only if its on the same strand as tr, hence this works. Also, no copy, since we do not modify f_align
    else:
        aligned=[tr['source_len']-pos for pos in reversed(tr['aligned'])]
        f_aligned=[[tr['source_len']-pos for pos in reversed(f['aligned'])] for f in fusions]

    used=set()# index of chimeric alignments explainable by splicing
    potential={i for i,f in enumerate(fusions) if f['chr']==chrom and f["strand"]==strand}#potentially explainable by splicing
    potential_pre = {i for i in potential if 0 < exons[0][0]-fusions[i]['exons'][-1][1] <= max_intron and aligned[0] > f_aligned[i][0]}#upstream primary transcript (wrt genome fwd strand) 
    potential_post= {i for i in potential if 0 < fusions[i]['exons'][0][0]-exons[-1][1] <= max_intron and aligned[1] < f_aligned[i][1]}#downstream primary transcript 
    found=True
    while found:
        found=False
        if potential_pre: #look upstream
            i=min(potential_pre , key=lambda idx: abs(aligned[0]-fusions[idx]['aligned'][1]) ) #todo: tiebreaker shortest intron length?
            delta=aligned[0]-f_aligned[i][1]
            if abs(delta) < tr_dist: #ideally this is 0
                intron_len=exons[0][0]-fusions[i]['exons'][-1][1]
                found=True
                #potential_pre.remove(i)
                used.add(i)
                cigar=merge_cigar(fusions[i]['cigar'],cigar, delta, intron_len) 
                aligned[0]=f_aligned[i][0] #update aligned
                idx=len(fusions[i]['exons'])-1
                exons=fusions[i]['exons']+exons #update exons
                if delta<0:
                    exons[idx][1]-=delta #overlapping alignments, extend one exon by abs(delta)
                potential_pre = {i for i in potential if 0 < exons[0][0]-fusions[i]['exons'][-1][1] <= max_intron and aligned[0] > f_aligned[i][0]}
    
        if potential_post and aligned[1]<tr['source_len']: #look downstream
            i=min(potential_post, key=lambda idx: abs(fusions[idx]['aligned'][0]-aligned[1]) ) #minimum gap (>0) or overlap (<0)
            delta=f_aligned[i][0]-aligned[1]
            if abs(delta) < tr_dist: #ideally this is 0
                intron_len=fusions[i]['exons'][0][0]-exons[-1][1]
                found=True
                #potential_post.remove(i)
                used.add(i)
                cigar=merge_cigar(cigar,fusions[i]['cigar'], delta, intron_len)
                aligned[1]=f_aligned[i][1]#update aligned  
                idx=len(exons)  -1           
                exons=exons+fusions[i]['exons']#update exons
                if delta<0:
                    exons[idx][1]-=delta #overlapping alignments 
                potential_post={i for i in potential if 0 < fusions[i]['exons'][0][0]-exons[-1][1] <= max_intron and aligned[1] < f_aligned[i][1]} #potentially explainable by splicing
    if used: #this is only needed for debugging, consider commenting this out. Its a bit inconsistant, since it modifies tr
        tr['chimeric_alignment_intron']=[(tr['cigar'],tr['aligned'], tr['exons'])]+[(fusions[f]['cigar'],fusions[f]['aligned'], fusions[f]['exons']) for f in used]
    if strand=='-':
        aligned=[tr['source_len']-pos for pos in reversed(aligned)]
    return exons,cigar,aligned, [f for i,f in enumerate(fusions) if i not in used]

def merge_cigar(cigar1, cigar2, delta,intron_len):
    'merge the cigar of two chimeric alignments which can be explained by splicing'
    idx1=next(i for i,c in enumerate(reversed(cigar1)) if c in 'M=X' ) #first index relevant cigar char from right
    for i,c in enumerate(cigar2):
        if not c.isdigit():
            if c in 'M=X':
                break
            idx2=i #last irrelevant (e.g. non reference) cigar char
    #this could lead to illegal cigar strings if there is a deletion directly after the soft split (which does not make sense, so i do not fix it here)
    if delta<0:
        gap=f'{-delta}D'
    elif delta>0:
        gap=f'{delta}I'
    else:
        gap=''
    return f'{cigar1[:-idx1]}{gap}{intron_len+(delta if delta<0 else 0)}N{cigar2[idx2+1:]}'
    


def get_transcript_sequence(chrom,exons,genome_fh):
    'construct the transcript sequence from genomic position and geneome sequence'
    #todo reverse complement if on (-) strand?
    seq=[genome_fh.fetch(chrom, *e) for e in exons]
    return ''.join(seq)


def get_mutations(cigartuples, seq, ref, chrom,ref_start):
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
            mutations.append((ref_pos,ref_base, alt_base))
        if cigar[0] in (0, 2, 3, 7, 8):  # MDN=X -> move forward on reference
            ref_pos += cigar[1]
        if cigar[0] in (0, 1, 4,7, 8):  # MIS=X -> move forward on seq
            seq_pos += cigar[1]
    return mutations        

def cigarstring2tuples(cigarstring):
    'converts cigar strings to cigartuples'
    #this is mainly for testing and debugging
    tuples=list()
    n=0
    for c in cigarstring:
        if c.isdigit():
            n*=10
            n+=ord(c)-48
        else:
            tuples.append((CIGAR_DICT[c],n))
            n=0
    return tuples

def find_orf(seq, min_len=3):
    'look for open reading in "seq" frames and return the longest'
    # todo: incorperate codon usage to filter false positives
    # assuming seq is upper case and stranded
    # this is quite inefficient, can be improved:
    #end=[0,0,0]
    #find start
    #if end[start%3]<start:
    #find next end[start%3]
    #if largest/best: save
    #repeat until end of seq -len(best)
    start_codon=[m.start() for m in re.finditer('ATG', seq)] #neglect rare alternative start codons: TTG and CTG
    stop_codon=sorted([m.start() for m in re.finditer('T[AA|AG|GA]', seq)])
    orf=list()
    for start in start_codon:
        try:
            stop=next(stop+3 for stop in stop_codon if stop>=start+min_len and stop%3 == start%3)
        except StopIteration:
            continue
        else:
            orf.append((start,stop))
    #orf.sort(key=lambda x:[x[0]-x[1]])
    if orf:
        return max(orf,key=lambda x:x[1]-x[0])
    return None




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
            # log.debug(line)
            try:
                gtf_id = info['transcript_id']
                _=exons.setdefault(gtf_id, list()).append((start, end))
            except KeyError:  # should not happen if GTF is OK
                warnings.warn("gtf format error: exon without transcript_id. Skipping line\n"+line)
        elif ls[2] == 'gene':
            if 'gene_id' not in info:
                warnings.warn("gtf format error: gene without gene_id. Skipping line\n"+line)
            else:
                info=_prepare_gene_info(info, ls[0], ls[6])
                new_gene=Gene(start, end, info, transcriptome)
                genes[ls[0]].add(new_gene)
                gene_dict[info['gene_id']]=new_gene
        elif ls[2] == 'transcript':
            try:
                _=transcripts.setdefault(info["gene_id"], list()).append(info["transcript_id"])
            except KeyError:
                warnings.warn("gtf format errror: transcript must have gene_id and transcript_id, skipping line\n"+line )
      
        else:
            skipped.add(ls[2])
    if len(skipped)>0:
        log.info('skipped the following categories: {}'.format(skipped))
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

def _set_alias(d,alias):
    for pref,alt in alias.items():
        alt=[a for a in alt if a in d]
        if pref not in d:
            try:
                d[pref]=next(d[a] for a in alt)
            except StopIteration:
                log.error(f'did not find alternative for {pref}- suggested terms are {alt}, but have only those keys:{list(d)}')
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
            warnings.warn("GFF format error in infos (should be ; seperated key=value pairs). Skipping line: "+line)
        start, end = [int(i) for i in ls[3:5]]
        start -= 1  # to make 0 based
        if ls[2] == "exon":
            try:
                gtf_id = info['Parent']
                exons.setdefault(gtf_id, list()).append((start, end))
            except KeyError:  # should not happen if GTF is OK
                warnings.warn("GFF format error: no parent found for exon. Skipping line: "+line)
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
        log.info('skipped the following categories: {}'.format(skipped))
    # sort the exons
    log.info('sorting exon positions...')
    for tid in tqdm(exons):
        exons[tid].sort()
    missed_genes=[gid for gid in transcripts.keys() if gid not in gene_set]
    if missed_genes:        
        #log.debug('/n'.join(gid+str(tr) for gid, tr in missed_genes.items()))
        notfound=len(missed_genes)
        found=sum((len(t) for t in genes.values()) )
        log.warning('Missing genes! Found gene information in categories {} for {}/{} genes'.format(gene_categories,found, found+notfound))
    log.info('building gene data structure...')
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

           

#utils

def filter_function(argnames, expression):
    'converts a string e.g. "all x[0]/x[1]>3" into a function'
    return eval (f'lambda {",".join(arg+"=None" for arg in argnames)}: bool({expression})\n',{},{})
    


def overlap(r1, r2):
    "check the overlap of two intervals"
    # assuming start < end
    if r1[1] < r2[0] or r2[1] < r1[0]:
        return False
    else:
        return True



def get_intersects(tr1, tr2):
    "get the number of intersecting splice sites and intersecting bases of two transcripts"
    tr1_enum = enumerate(tr1)
    try:
        j,tr1_exon = next(tr1_enum)
    except StopIteration:
        return 0, 0
    sjintersect = 0
    intersect = 0
    for i, tr2_exon in enumerate(tr2):
        while tr1_exon[0] < tr2_exon[1]:
            if tr2_exon[0] == tr1_exon[0] and i > 0 and j>0:  # neglegt TSS and polyA
                sjintersect += 1
            if tr2_exon[1] == tr1_exon[1] and i < len(tr2)-1 and j <len(tr1)-1:
                sjintersect += 1
            if overlap(tr1_exon, tr2_exon):
                # region intersect
                i_end = min(tr1_exon[1], tr2_exon[1])
                i_start = max(tr1_exon[0], tr2_exon[0])
                intersect += (i_end-i_start)
            try:
                j,tr1_exon = next(tr1_enum)
            except StopIteration:  # tr1 is at end
                return sjintersect, intersect
    # tr2 is at end
    return sjintersect, intersect

def get_relation(exons1, exons2):
    ''' returns a list of length of exons1.
     each element represents a exon from set 1 and contains a list of pairs, denoting corresponding exons from set 2
     the first element of the pair is the index of the corresponding set 2 exon,
     the second element is in [0,1,2,3] and encodes the splice correspondence: 0 overlap only, 1: same splice donor, 2: same splice acceptor, 3: identical'''
    relation = [[] for _ in exons1]
    enum2 = enumerate(exons2)
    i, e2 = next(enum2)
    for j, e1 in enumerate(exons1):
        while e1[1] > e2[0]:
            ol = overlap(e1, e2)-1
            if ol >= 0:
                if e1[1] == e2[1]:
                    ol += 1
                if e1[0] == e2[0]:
                    ol += 2
                relation[j].append((i, ol))
            if e2[1]<=e1[1]:
                try:
                    i, e2 = next(enum2)
                except StopIteration:  # exons2 is at end
                    return relation
            else: break
    return relation


def is_truncation(long_tr, short_tr): 
    'checks whether short a truncated version of long is'
    relation=get_relation(short_tr, long_tr)
    if any(len(r) != 1 for r in relation):
        return False ,0,0
    if any(r[0][1]!=3 for r in relation[1:-1]):
        return False,0,0
    if len(relation)>0:
        if relation[0][0][1] & 1 ==0: #only for multi exon transcripts
            return False,0,0
        if relation[-1][0][1] & 2 ==0:
            return False,0,0
    if relation[0][0][0] > 0 and short_tr[0][0]<long_tr[relation[0][0][0]][0]:
        return False,0,0
    if relation[-1][0][0] < len(long_tr)-1 and short_tr[-1][1]>long_tr[relation[-1][0][0]][1]:
        return False,0,0    
    if any(relation[i][0][0]-relation[i+1][0][0]!=-1 for i in range(len(relation)-1)):
        return False, 0,0
    delta1=sum( [long_tr[i][1] -long_tr[i][0] for i in range(relation[0][0][0])]) #exons only in long
    delta1+=short_tr[0][0]-long_tr[relation[0][0][0]][0] #first overlapping exon
    delta2=sum( [long_tr[i][1] -long_tr[i][0] for i in range(relation[-1][0][0]+1, len(long_tr))])#exons only in long
    delta2+=long_tr[relation[-1][0][0]][1]-short_tr[-1][1]#last overlapping exon
    return True, delta1, delta2

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

def junctions_from_cigar(cigartuples, offset):
    'returns the exon positions'
    exons = list([[offset, offset]])
    for cigar in cigartuples:
        if cigar[0] == 3:  # N ->  Splice junction
            pos = exons[-1][1]+cigar[1]
            if exons[-1][0]==exons[-1][1]: 
                # delete zero length exons 
                # (may occur if insertion within intron, e.g. 10M100N10I100N10M)
                del exons[-1]
            exons.append([pos, pos])
        elif cigar[0] in (0, 2, 7, 8):  # MD=X -> move forward on reference
            exons[-1][1] += cigar[1]
    if exons[-1][0]==exons[-1][1]: #delete 0 length exons at the end
        del exons[-1]
    return exons

def is_same_transcript(tr1, tr2):
    'checks whether tr1 and tr2 are the same gene by calculating intersection over union of the intersects'

def is_same_gene(tr1, tr2, spj_iou_th=0, reg_iou_th=.5):
    'checks whether tr1 and tr2 are the same gene by calculating intersection over union of the intersects'
    #current default definition of "same gene": at least one shared splice site
    # or more than 50% exonic overlap
    spj_i, reg_i=get_intersects(tr1, tr2)
    total_spj=(len(tr1)+len(tr2)-2)*2
    spj_iou=spj_i/(total_spj-spj_i) if total_spj>0 else 0
    if spj_iou>spj_iou_th:
        return True
    total_len=sum([e[1]-e[0] for e in tr2+tr1])
    reg_iou=reg_i/(total_len-reg_i)
    if reg_iou>reg_iou_th:
        return True
    return False

def splice_identical(tr1, tr2):
    #all splice sites are equal
    if len(tr1) != len(tr2):  # different number of exons
        return False
    if len(tr1) == 1 and overlap(tr1[0], tr2[0]):  # single exon genes
        return True
    if tr1[0][1] != tr2[0][1] or tr1[-1][0] != tr2[-1][0]:  # check first and last exons
        return False
    for e1, e2 in zip(tr1[1:-1], tr2[1:-1]):  # check other exons
        if e1[0] != e2[0] or e1[1] != e2[1]:
            return False
    return True
