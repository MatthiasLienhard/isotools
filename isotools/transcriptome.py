#!/usr/bin/python3
import os
import os.path

import itertools
import warnings
from itertools import combinations, product
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
        'CLIPPED_ALIGNMENT':'clipped>50',#more then 50 bases or 20% clipped
        'A_CONTENT':'downstream_A_content>.5', #more than 50% a
        'RTTS':'template_switching and any(ts[2]<=10 and ts[2]/(ts[2]+ts[3])<.2 for ts in template_switching)', #less than 10 reads and less then 20% of total reads for at least one junction
        'NONCANONICAL_SPLICING':'noncanonical_splicing',
        'NOVEL_TRANSCRIPT':'annotation is None or "splice_identical" not in annotation["as"]',
        'TRUNCATION':'truncated',
        'REFERENCE':'annotation',# and "splice_identical" in annotation',
        'UNSPLICED':'len(exons)==1','MULTIEXON':'len(exons)>1'}
        
class Transcriptome:
    def __init__(self, file_name, file_format='auto',**kwargs):        
        self.data,self.infos=import_transcripts(file_name, file_format=file_format, **kwargs)
        if 'file_name' not in self.infos:
            self.infos['file_name']=file_name
        self.make_index()
        
    def save(self, fn=None):
        if fn is None:
            fn=self.infos['file_name']+'.isotools.pkl'
        log.info('saving transcriptome to '+fn)
        pickle.dump((self.data,self.infos), open(fn, 'wb'))
    
    def load(self, fn=None):
        if fn is None:
            fn=self.infos['file_name']+'.isotools.pkl'
        log.info('loading transcriptome from '+fn)
        self.data, self.infos=pickle.load(open(fn, 'rb'))
        
    def write_gtf(self, fn, source='isotools',use_gene_name=False,  include=None, remove=None):     
        with open(fn, 'w') as f:     
            for gene in tqdm(self):
                lines=gene.to_gtf(source=source,use_gene_name=use_gene_name, include=include, remove=remove)
                if lines:
                    _=f.write('\n'.join( ('\t'.join(str(field) for field in line) for line in lines) )+'\n')



    def make_index(self):
        idx=dict()
        for g in self:
            idx[g.name] = g
            idx[g.id]=g
        self._idx=idx
        
    def __getitem__(self, key):
        return self._idx[key]

    def __len__(self):
        return self.n_genes
    
    def remove_chromosome(self, chromosome):
        for n in (g.id for g in self.data[chromosome]):
            self._idx[n]={g for g in self._idx[n] if g.chrom != chromosome}
            if not self._idx[n]:
                del self._idx[n]
        del self.data[chromosome]
    
    def add_illumina_coverage(self, illumina_fn, load=False):
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
        #populates transcript['biases'] which is the bases for the filter
        with FastaFile(genome_fn) as genome_fh:
            for g in tqdm(self):
                if 'splice_graph' not in g.data:
                    g.data['splice_graph']=isotools.splice_graph.SpliceGraph(g)
                ts_candidates=g.data['splice_graph'].find_ts_candidates()
                for start, end, js, ls, idx in ts_candidates:
                    for tr in (g.transcripts[i] for i in idx):
                        tr.setdefault('template_switching',[]).append((start, end, js, ls))
                g.add_direct_repeat_len(genome_fh) 
                g.add_noncanonical_splicing(genome_fh)
                g.add_threeprime_a_content(genome_fh)
                g.add_truncations()
        self.infos['biases']=True # flag to check that the function was called

    def add_filter(self, transcript_filter={},gene_filter={}):
        #possible filter flags: 'A_CONTENT','RTTS','NONCANONICAL_SPLICING','NOVEL_GENE','NOVEL_TRANSCRIPT','TRUNCATION'
        # a_th=.5, rtts_maxcov=10, rtts_ratio=5, min_alignment_cov=.1
        #biases are evaluated with respect to specific thresholds
        gene_attributes={k for g in self for k in g.data.keys() }
        tr_attributes={k for g in self for tr in g.transcripts.values() for k in tr.keys() }
        tr_attributes.add('filter')
        if not gene_filter and not transcript_filter:
            gene_filter=default_gene_filter
            transcript_filter=default_transcript_filter
        gene_ffun={label:filter_function(gene_attributes, fun) for label,fun in gene_filter.items()}
        tr_ffun={label:filter_function(tr_attributes, fun) for label,fun in transcript_filter.items()}
        for g in tqdm(self):
            g.add_filter(tr_ffun,gene_ffun)
        self.infos['filter']={'gene_filter':gene_filter, 'transcript_filter':transcript_filter}

    @property 
    def runs(self):
        return dict(zip(self.infos['sample_table'].name, self.infos['sample_table'].run))
        
    def get_sample_idx(self, group_column='name'):
        #returns a dict with group names as keys and index lists as values
        return self.infos['sample_table'].groupby(group_column).groups
        
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
    def chromosomes(self):
        return list(self.data)            

    def __str__(self):
        return '{} object with {} genes and {} transcripts'.format(type(self).__name__, self.n_genes, self.n_transcripts)
    
    def __repr__(self):
        return object.__repr__(self)

    def __iter__(self):
        return (gene for tree in self.data.values() for gene in tree)

    

    def annotate(self, reference,fuzzy_junction=5,novel_params=None):  
        genes=dict()
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
                        exons=tr['exons']
                        if chrom in reference.data:
                            ref_genes_ol = [rg for rg in reference.data[chrom][exons[0][0]: exons[-1][1]] if rg.strand==g.strand]
                            if len(ref_genes_ol)>0:                                
                                intersect=[rg.splice_graph.get_intersects(exons) for rg in ref_genes_ol]
                                best_idx ,(sj_i, base_i)= max(enumerate(intersect), key=lambda x:x[1])           #first criterion: sj overlap, second: base overlap                     
                                if base_i > 0 :    # at least some overlap
                                    ref_g = ref_genes_ol[best_idx]
                                    g.correct_fuzzy_junctions(trid,ref_g, fuzzy_junction)
                                    altsplice= ref_g.splice_graph.get_alternative_splicing(exons, g.strand)
                                    tr['annotation']={'ref_gene_name': ref_g.name, 'ref_gene_id':ref_g.id, 'sj_i': sj_i, 'base_i':base_i,'as':altsplice}
                                    #todo: what about fusion?  
                                    gene_idx.setdefault(ref_g.name,{'chr':chrom, 'ID':ref_g.id, 'Name': ref_g.name, 'strand':g.strand ,'transcripts': dict()})
                                    gene_idx[ref_g.name]['transcripts'][trid]=tr
                                    continue
                        tr['annotation']=None
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
                    genes[chrom].add(Gene(start, end, g))
        self.data=genes
        self.infos['annotation']=reference.infos['file_name']
        

    def iter_genes(self, region):
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
            yield g
    
    def iter_transcripts(self,region=None,include=None, remove=None):
        for g in self.iter_genes(region):
            for trid in g.filter_transcripts(include, remove):
                yield g,trid,g.transcripts[trid]


    def gene_table(self, region=None ): #ideas: filter, extra_columns
        colnames=['chr', 'start', 'end', 'strand', 'gene_name', 'n_transcripts']        
        rows=[(g.chrom, g.start, g.end, g.strand, g.id, g.n_transcripts) for g in  self.iter_genes(region)]
        df = pd.DataFrame(rows, columns=colnames)
        return(df)

    def transcript_table(self, region=None, extra_columns=None,  include=None, remove=None): 
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
    
class Coverage: #stores the illumina read coverage of a gene 
    #plan: make a binned version, or use run length encoding
    def __init__(self, cov, junctions, offset, chrom=None):
        self._cov=cov
        self._junctions=junctions
        self.reg=None if cov is None else (chrom,offset, offset+len(cov)) 
        self.bam_fn=None
    
    @classmethod
    def from_bam(cls, bam_fn,g, load=False):
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
        cov,junctions=cls._import_coverage(align_fh,(g.chrom, g.start,g.end))
        obj = cls.__new__(cls)  
        obj.__init__(cov,junctions, g.start )
        return obj

    @classmethod #this is slow - called only if coverage is requested
    def _import_coverage(cls,align_fh, reg):
        delta=np.zeros(reg[2]-reg[1])        
        junctions={}
        for read in align_fh.fetch(*reg):        
            exons=junctions_from_cigar(read.cigartuples,read.reference_start) #alternative: read.get_blocks() should be more efficient.. -todo: is it different?
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
        cov=np.cumsum(delta)   
        return cov,junctions
    
    def load(self):
        with AlignmentFile(self.bam_fn, 'rb') as align:
            log.info(f'Illumina coverage of region {self.reg[0]}:{self.reg[1]}-{self.reg[2]} is loaded from {self.bam_fn}') #info or debug?
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
        return self._cov

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
    def __str__(self):
        return('Gene {} {}({}), {} transcripts'.format(self.name, self.region, self.strand, self.n_transcripts))
    
    def __repr__(self):
        return object.__repr__(self)

    def add_filter(self, transcript_filter,gene_filter):   
        gene_tags=[name for name,fun in gene_filter.items() if  fun(**self.data)]
        for tr in self.transcripts.values():
            tr['filter']=gene_tags.copy()
            tr['filter'].extend([name for name,fun in transcript_filter.items() if fun(**tr)])
            if not tr['filter']:
                tr['filter']=['PASS']


    def filter_transcripts(self, include=None, remove=None):
        #include: transcripts must have at least one of the flags
        #remove: transcripts must not have one of the flags (priority)
        if include is None and remove is None: # return all
            for tr in self.transcripts.keys():
                yield tr
        else:
            for trid in self.transcripts.keys():            
                try:
                    select =  any(f in self.transcripts[trid]['filter'] for f in include)
                except TypeError:
                    select=True #None means include all
                try:
                    select &= all(f not in self.transcripts[trid]['filter'] for f in remove)
                except TypeError:
                    pass
                if select:
                    yield trid
        
    def correct_fuzzy_junctions(self,trid,ref_g, size):
        exons=self.transcripts[trid]['exons']
        shifts=ref_g.splice_graph.fuzzy_junction(exons,size)
        if shifts:
            if 'splice_graph' in self.data:
                del self.data['splice_graph']
            self.transcripts[trid]['fuzzy_junction']=shifts #save info for analysis of the effect
            for i,sh in shifts.items():
                exons[i][1]+=sh
                exons[i+1][0]+=sh

    def to_gtf(self, source='isoseq', use_gene_name=False, include=None, remove=None):
        
        info={'gene_id':self.name if use_gene_name else self.id}
        gene=(self.chrom, source, 'gene', self.start+1, self.end, '.',self.strand, '.', '; '.join('{} "{}"'.format(k,v) for k,v in info.items() ))
        lines=list()
        for trid in self.filter_transcripts(include, remove):
            tr=self.transcripts[trid]
            info['transcript_id']=trid
            #todo: print only relevant infos
            lines.append((self.chrom, source, 'transcript', tr['exons'][0][0]+1, tr['exons'][-1][1], '.',self.strand, '.', '; '.join('{} "{}"'.format(k,v) for k,v in info.items() if k != 'exon_id')))
            for enr, pos in enumerate(tr['exons']):
                info['exon_id']='{}_{}'.format(trid, enr)
                lines.append((self.chrom, source, 'exon', pos[0]+1, pos[1], '.',self.strand, '.', '; '.join('{} "{}"'.format(k,v) for k,v in info.items())))
        if lines:
            lines.append(gene)
        return(lines)
    


    def add_noncanonical_splicing(self, genome_fh):
        #for all transcripts, add the the index and sequence of noncanonical (i.e. not GT-AG)
        #i.e. save the dinucleotides of donor and aceptor as XX-YY string in the transcript biases dicts
        #noncanonical splicing might indicate technical artifacts (template switching, missalignment, ...)
        
        for tr in self.transcripts.values():
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
        #perform a global alignment of the sequences at splice sites and report the score. Direct repeats indicate template switching
        for tr in self.transcripts.values():
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
        #add the information of the genomic A content downstream the transcript. High values indicate genomic origin of the pacbio read        
        for tr in self.transcripts.values():
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
        #check for truncations and save indices of truncated genes in transcripts
        #dist5=-1, dist3=10
        is_truncated=set()
        for trid1, trid2 in combinations(self.transcripts.keys(), 2):
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
        #returns length of 5'UTR, coding sequence and 3'UTR
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
            raise
            
    @property
    def name(self):
        try:
            return self.data['Name']    
        except KeyError:
            pass
        try:
            return self.data['gene_name']    
        except KeyError:
            pass        
        
        return self.id

    #@id.setter
    #def id(self, new_id):
    #    self.data['gene_id']=new_id

    @property
    def strand(self):
        try:
            return self.data['strand']
        except KeyError:
            return {}

    @property
    def transcripts(self):
        try:
            return self.data['transcripts']
        except KeyError:
            return {}

    @property
    def n_transcripts(self):
        try:
            return len(self.data['transcripts'])
        except KeyError:
            return 0

    @property
    def splice_graph(self):
        try:
            return self.data['splice_graph']
        except KeyError:
            self.data['splice_graph']=isotools.splice_graph.SpliceGraph(self)
            return self.data['splice_graph']

    def __copy__(self):
        return Gene(self.start, self.end, self.data)        
        

    def __deepcopy__(self, memo):
        return Gene(self.start, self.end, copy.deepcopy(self.data, memo))
        

    def __reduce__(self):
        return Gene, (self.start, self.end, self.data)  

    def copy(self):
        return self.__copy__()

def genomic_position(exons, tr_pos):
    # gets a sorted list of transcript positions and finds the corresponding genomic positions
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
                


def import_transcripts(fn, file_format='auto', **kwargs):
    if file_format=='auto':        
        file_format=os.path.splitext(fn)[1].lstrip('.')
        if file_format=='gz':
            file_format=os.path.splitext(fn[:-3])[1].lstrip('.')
    log.info(f'importing transcriptome from {file_format} file {fn}')
    if file_format == 'gtf':
        genes,infos= import_gtf_transcripts(fn,  **kwargs)
    elif file_format in ('gff', 'gff3'):
        genes,infos= import_gff_transcripts(fn,  **kwargs)
    elif file_format == 'bam':
        genes,infos= import_pacbio_transcripts(fn,  **kwargs)
    elif file_format == 'pkl':
        genes,infos= pickle.load(open(fn, 'rb'))
    else:
        raise ValueError('unsupportet file format: {}'.format(file_format))    
    return genes,infos

def import_pacbio_transcripts(fn,sample_table,genome_fn=None, chromosomes=None):
        n_secondary=0
        neglected_runs=set()
        #runs=pd.DataFrame(runs)
        assert isinstance(sample_table, pd.DataFrame) and 'name' in sample_table.columns and 'run' in sample_table.columns, '\"sample_table\" must be a pandas DataFrame with columns \"name\" and \"run\"'
        sample_table.reset_index(drop=True)# make sure the index is 0..n
        run_idx={v:k for k,v in sample_table.run.items()}
        #with FastaFile(genome_fn) if genome_fn is not None else None as genome_fh: #none has no open, does not work
        try:
            genome_fh=FastaFile(genome_fn) if genome_fn is not None else None
            with AlignmentFile(fn, "rb") as align:        
                stats = align.get_index_statistics()
                # try catch if sam/ no index /not pacbio?
                total_reads = sum([s.mapped for s in stats])
                if chromosomes is None:
                    chromosomes=align.references
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
                    cds=find_orf(read.query_sequence)
                    
                    #transcript infos
                    tr={'exons':exons, 'source_len':len(read.query_sequence), 'cigar':read.cigarstring, 
                        'clipped': sum(n for i,n in read.cigartuples if i in (4,5)),
                        'coverage':cov, 'CDS':cds} 
                    if genome_fh is not None:
                        tr['mutations']=get_mutations(read.cigartuples, read.query_sequence, genome_fh, chrom,read.reference_start)
                    #gene infos
                    info={'strand':strand, 'chr':chrom,
                            'ID':read.query_name,'transcripts':{read.query_name:tr}}
                    
                    if 'SA' in tags: #part of a fusion gene   
                        if read.flag & 0x800: #secondary part
                            fusion.setdefault(read.query_name, []).append({'chr':chrom,'exons':tr['exons'],'cigar':tr['cigar'],'strand':strand}) #cigar is required to sort out the order 
                        else: #primary part
                            fusion_primary[read.query_name]=tr
                    if not read.flag & 0x800:
                        genes[chrom].add(Gene(exons[0][0],exons[-1][1],info))
        except:
            if genome_fh is not None:
                genome_fh.close()
            raise
        if genome_fh is not None:
            genome_fh.close()
        if neglected_runs:
            log.warning(f'found the follwing runs in pacbio file, which where not specified in the run table: {neglected_runs}')    
        #link secondary alignments (chimeric genes)
        for tid,primary in fusion_primary.items():
            primary['fusion']=fusion[tid]
        if n_secondary>0:
            log.info(f"skipped {n_secondary} secondary transcripts {n_secondary/total_reads*100}%")
        return genes,{'sample_table':sample_table}


def get_transcript_sequence(chrom,exons,genome_fh):
    seq=[genome_fh.fetch(chrom, *e) for e in exons]
    return ''.join(seq)


def get_mutations(cigartuples, seq, ref, chrom,ref_start):
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

def find_orf(seq, min_len=3):
    # look for open reading frames and return the longest
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




def import_gtf_transcripts(fn, genes=None,chromosomes=None):
    gtf = TabixFile(fn)   
    exons = dict()  # transcript id -> exons
    transcripts = dict()  # gene_id -> transcripts
    skipped = set()
    if genes is None:
        genes=dict()  
    elif(isinstance(genes, Transcriptome)) :
        genes=genes.data
    gene_dict={g.data['name']:g for t in genes.values() for g in t}
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
                info=prepare_gene_info(info, ls[0], ls[6])
                new_gene=Gene(start, end, info)
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
            gene.data.setdefault('transcripts', {})
            for t_id in t_ids:
                try:
                    gene.data['transcripts'][t_id] = {'exons':exons[t_id]}
                except KeyError:
                    # genes without transcripts get a single exons transcript
                    gene.data['transcripts'] = {t_id: {'exons':[tuple(gene[:2])]}}                
    return genes,dict()

def prepare_gene_info(info,chrom, strand, modify=True):
    if not modify:
        info=copy.deepcopy(info)
    if 'name' not in info:
        if 'Name' in info:
            info['name']=info['gene_id']    
            del info['Name']
        else:
            info['name']=info['gene_id']
            del info['gene_id']
    info['strand'] =strand
    info['chr'] = chrom
    return info

def import_gff_transcripts(fn, chromosomes=None, gene_categories=['gene']):
    # returns a dict interval trees for the genes, each containing the splice graph
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
            gene_set.add(info['ID'])
            genes[chrom].add(Gene(start,end,info))
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
        log.warning('found specific gene information with category {} for {}/{} genes'.format(gene_categories,found, found+notfound))
    log.info('building gene data structure...')
    # add transcripts to genes
    for chrom in genes:
        for gene in genes[chrom]:
            g_id = gene.data['ID']
            tr=transcripts.get(g_id,  {g_id:{}})
            for t_id,tr_info in tr.items():
                try:
                    tr_info['exons']=exons[t_id]
                except KeyError:
                    # genes without exons get a single exons transcript
                    tr_info['exons']=[tuple(gene[:2])]
                #add cds
                if t_id in cds_start and t_id in cds_stop:
                    tr_info['CDS']=(cds_start[t_id], cds_stop[t_id]) if cds_start[t_id]< cds_stop[t_id] else (cds_stop[t_id],cds_start[t_id])
                gene.data.setdefault('transcripts', dict())[t_id]=tr_info
    return genes,dict()

def get_gff_chrom_dict(gff, chromosomes):
    # fetch chromosome ids - in case they use ids in gff for the chormosomes
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

def collapse_transcripts_to_genes( transcripts, spj_iou_th=0, reg_iou_th=.5, gene_prefix='PB_novel_', offset=0):
    #transcripts is an interval tree of genes (e.g. from one chromosome)
    n_novel=0
    idx=dict()
    merge=dict()
    for g in transcripts:
        for trid,tr in g.transcripts.items():
            idx[trid]=tr
            merge_trid={trid}
            candidates=[c for c in transcripts.overlap(g.start, g.end) if c.strand==g.strand]
            for ol_trid,ol_tr in [(k,v) for c in candidates for k,v in c.transcripts.items()]:
                if ol_trid in merge_trid  or ol_trid not in merge: # proceed if already added or not seen in the main loop so far
                    continue
                if is_same_gene(tr['exons'], ol_tr['exons'],spj_iou_th, reg_iou_th):
                    #add all transcripts of overlapping
                    merge_trid.update(merge[ol_trid])
            #update all overlapping (add the current to them)
            for ol_trid in merge_trid:
                merge[ol_trid]=merge_trid
    genes=IntervalTree()       
    seen=set()
    for g in transcripts:
        for trid,tr in g.transcripts.items():
            if trid in seen:
                continue
            n_novel+=1
            new_transcripts={tid:idx[tid] for tid in merge[trid]}
            seen.update(merge[trid])
            start=min(tr['exons'][0][0] for tr in new_transcripts.values())
            end=max(tr['exons'][-1][1] for tr in new_transcripts.values())
            new_data={'chr':g.chrom, 'ID':f'{gene_prefix}{n_novel+offset}', 'strand':g.strand, 'transcripts':new_transcripts}
            genes.add(Gene(start,end,new_data ))
    return genes
               

#utils

def filter_function(argnames, expression):
    #converts a string e.g. 'all x[0]/x[1]>3' into a function
    #todo: is it save??
    return eval (f'lambda {",".join(arg+"=None" for arg in argnames)}: bool({expression})\n',{},{})
    


def overlap(r1, r2):
    # do the two intervals overlap
    # assuming start < end
    if r1[1] < r2[0] or r2[1] < r1[0]:
        return False
    else:
        return True



def get_intersects(tr1, tr2):
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
    # returns a list of length of exons1.
    # each element represents a exon from set 1 and contains a list of pairs, denoting corresponding exons from set 2
    # the first element of the pair is the index of the corresponding set 2 exon,
    # the second element is in [0,1,2,3] and encodes the splice correspondence: 0 overlap only, 1: same splice donor, 2: same splice acceptor, 3: identical
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
    #is short a truncated version of long?
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

def junctions_from_cigar(cigartuples, offset):
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
    return exons

def is_same_gene(tr1, tr2, spj_iou_th=0, reg_iou_th=.5):
    # current default definition of "same gene": at least one shared splice site
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
