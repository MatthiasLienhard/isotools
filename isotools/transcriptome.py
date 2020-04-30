#!/usr/bin/python3
import os
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
import logging
import copy
import pickle
import operator as o

log=logging.getLogger(__name__)
log.setLevel(logging.INFO)
log_format=logging.Formatter('%(levelname)s: [%(asctime)s] %(name)s:%(message)s')
#log_file=logging.FileHandler('logfile.txt')
log_stream=logging.StreamHandler()
log_stream.setFormatter(log_format)
log.addHandler(log_stream)
#a_th=.5, rtts_maxcov=10, rtts_ratio=5, min_alignment_cov=.1
default_gene_filter={'NOVEL_GENE':'all(tr["annotation"] for tr in transcripts)'}
default_transcript_filter={
        'CLIPPED_ALIGNMENT':'clipped>50',#more then 50 bases or 20% clipped
        'A_CONTENT':'downstream_A_content>.5', #more than 50% a
        'RTTS':'any(ts[2]<=10 and ts[2]/(ts[2]+ts[3])<.2 for ts in template_switching)', #less than 10 reads and less then 20% of total reads for at least one junction
        'NONCANONICAL_SPLICING':'noncanonical_splicing',
        'NOVEL_TRANSCRIPT':'not annotation or "splice_identical" not in annotation["sType"]',
        'TRUNCATION':'truncated'}
        
class Transcriptome:
    def __init__(self, file_name, chromosomes=None, file_format='auto', groups=None):        
        self.data,self.infos=import_transcripts(file_name, chromosomes=chromosomes, file_format=file_format)
        if groups is not None:
            self.groups=groups
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

    '''
    def write_table(self, fn, source='isotools'):     
        header=['id','name', 'chrom','strand','txStart','txEnd','exonCount','exonStarts','exonEnds']
        with open(fn, 'w') as f:       
            print('\t'.join(header))
            for gene in tqdm(self):
                lines=gene.to_table()
                print('\n'.join( ('\t'.join(str(field) for field in line) for line in lines) ),file=f)
    '''

    def make_index(self,use_name=False):
        idx=dict()
        for g in self:
            idx[g.name if use_name else g.id]=g
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
        return(self.infos['runs'])
        
    @property 
    def groups(self):
        try:
            return(self.infos['groups'])
        except KeyError:
            return({r:[r] for r in self.infos['runs']})

    @groups.setter
    def groups(self, grp):
        assert all(r in self.infos['runs'] for g in grp for r in g)
        self.infos['groups']=grp

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

    

    def annotate(self, reference,novel_params=None):  
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
                        tr['annotation']=get_support(tr['exons'], reference, chrom=chrom,is_reverse=g.strand == '-')
                        #todo: what about fusion?                
                        if tr['annotation'] is not None:
                            gname=tr['annotation']['ref_gene_name']
                            gid=tr['annotation']['ref_gene_id']          
                            if gname not in gene_idx:
                                gene_idx[gname]={'chr':chrom, 'ID':gid, 'Name': gname, 'strand':g.strand ,'transcripts': dict()}
                            gene_idx[gname]['transcripts'][trid]=tr
                        else:                         
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
            #todo:merge and add novel genes
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
        extracolnames={'alt_splice':('splice_IoU', 'base_IoU','splice_type')}
        if 'runs' in self.infos:
            extracolnames['coverage']=(f'coverage_{r}' for r in self.infos['runs'])
        colnames=['chr', 'gene_start', 'gene_end','strand', 'gene_id','gene_name' ,'transcript_name']     + [n for w in extra_columns for n in (extracolnames[w] if w in extracolnames else (w,))]
        rows=[]
        for g,trid, tr in self.iter_transcripts(region, include,remove):                
                rows.append((g.chrom, g.start, g.end, g.strand,g.id, g.name, trid)+g.get_values(trid,extra_columns))
        df = pd.DataFrame(rows, columns=colnames)
        return(df)
    
    
        
#def default_false(fun,**args):
#    try:
#        return fun(**args)
#    except NameError:
#        return False

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

    '''
    def add_filter(self, a_th=.5, rtts_maxcov=10, rtts_ratio=5, min_alignment_cov=.1):   
       novel_gene=all(tr['annotation'] is None for tr in self.transcripts.values())
        for tr in self.transcripts.values():
            flags=[]
            if abs(tr['source_len'] - sum(e-s for s,e in tr['exons']))/tr['source_len']>min_alignment_cov:
                flags.append('CLIPPED_ALIGNMENT')
            if tr['biases']['downstream_A_content']>a_th:
                flags.append('A_CONTENT')                
            if 'template_switching' in tr['biases']:
                for ts in tr['biases']['template_switching']:
                    s, l=(int(i) for i in ts[ts.rfind(':')+1:].split('/',1))
                    if s<rtts_maxcov and l/s>rtts_ratio:
                        flags.append('RTTS')
                        break
            if 'noncanonical_splicing' in tr['biases']:
                flags.append('NONCANONICAL_SPLICING')
            if novel_gene:
                flags.append('NOVEL_GENE')
            elif tr['annotation'] is None or 'splice_identical' not in tr['annotation']['sType']:
                flags.append('NOVEL_TRANSCRIPT')
            if 'truncated' in tr['biases']:    
                flags.append('TRUNCATION')
            tr['filter']=flags
    '''        

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
    
    ''' does not work so far (depends on unimplemented splice graph functionality)
    def unify_junctions(self, fuzzy_junction=5, reference=None):
        if 'splice_graph' not in self.data:
            self.data['splice_graph']=isotools.splice_graph.SpliceGraph((tr['exons'] for tr in self.transcripts), self.end)
        new_start, new_end=self.data['splice_graph'].unify_junctions(fuzzy_junction)  
        for tr in self.transcripts:
            try:
                tr['exons']=[(new_start[e[0]], new_end[e[1]]) for e in tr['exons']] #todo: simply use zip?
            except KeyError:
                logging.info('no new values for transcript {}: {}\nnew starts:{}\nnew ends:{}'.format(tr['name'],tr['exons'],new_start, new_end))
                raise
    '''

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
        elif what=='alt_splice':
            sel=['sjIoU','exIoU', 'sType']
            support=self.transcripts[trid]['annotation']
            if support is None:
                return ('NA',)*len(sel)
            else:
                #vals=support[n] if n in support else 'NA' for n in sel
                stype=support['sType']
                if isinstance(stype, dict):
                    type_string=';'.join('{}:{}'.format(k,v) for k,v in stype.items())
                else:
                    type_string=';'.join(str(x) for x in stype)
                vals=(support['sjIoU'],support['exIoU'],type_string)
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
    def start(self):
        return self.begin
    
    @property
    def region(self):
        try:
            return '{}:{}-{}'.format(self.chrom,self.start, self.end )
        except KeyError:
            raise


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
            return 'NA'

    @property
    def transcripts(self):
        try:
            return self.data['transcripts']
        except KeyError:
            return 'NA'

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

'''
def merge_genes(gene, candidates,spj_iou_th, reg_iou_th):
    same=set()
    new_gene=False
    for ol_gene in candidates:
        if ol_gene.strand != gene.strand or gene == ol_gene :
            continue
        if any(is_same_gene(tr1['exons'], tr2['exons'],spj_iou_th, reg_iou_th) for tr1, tr2 in itertools.product(gene.transcripts(), ol_gene.transcripts())):
            same.add(ol_gene)
    if len(same)>0:
        same.add(gene)
        new_data={'strand': gene.strand, 'chr': gene.chrom, 'source':{g.id for g in same}, 'transcripts':dict(ChainMap(*(g.transcripts for g in same)))}
        #todo merge more properties
        new_gene=Gene(min(g.start for g in same), max(g.end for g in same),new_data )
        assert len(new_gene.transcripts)==sum(len(g.transcripts) for g in same), 'new gene has different number of transcripts compared to source!\nsource:\n{}\nmerged\n{}'.format(same, new_gene)
    return new_gene, same

def merge_transcripts(tr1, tr2):
    #print('merging {} {}'.format(tr1, tr2))
    merged={'exons':list()}
    e2iter=iter(tr2['exons'])
    e1enum=enumerate(tr1['exons'])
    e2=next(e2iter)
    for i,e1 in e1enum:
        if overlap(e1,e2):
            merged['exons'].append( ( min(e1[0], e2[0]), max(e1[1], e2[1]) ) )
        elif e1[0]<e2[0]:
            merged['exons'].append(e1)
            continue
        else:
            merged['exons'].append(e2)
        try:
            e2=next(e2iter)
        except StopIteration:
            merged['exons']+=[e for j,e in e1enum if j> i]
            break
    else:
        merged['exons'].append(e2)
        merged['exons']+=list(e2iter)

    for n in tr1.keys()|tr2.keys():
        if n not in ['exons', 'annotation']:
            merged[n]=tr1.get(n,[])+tr2.get(n,[]) #carefull, might not be possible for all datatypes

    return merged
'''        




#io
def get_read_sequence(bam_fn,reg=None, name=None):
    with AlignmentFile(bam_fn, "rb") as align:
        for read in align.fetch(region=reg):
            chrom = read.reference_name
            strand = '-' if read.is_reverse else '+'
            exons = junctions_from_cigar(read.cigartuples,read.reference_start)
            if name is None or name == read.query_name:
                yield  read.query_name, read.query_sequence, chrom,strand, exons



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

def import_pacbio_transcripts(fn,chromosomes=None):
        n_secondary=0
        with AlignmentFile(fn, "rb") as align:        
            stats = align.get_index_statistics()
            # try catch if sam/ no index /not pacbio?
            total_reads = sum([s.mapped for s in stats])
            if chromosomes is None:
                chromosomes=align.references
            runs=set()
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
                tr={'exons':exons, 'source_len':len(read.seq), 'cigar':read.cigarstring, 'clipped': sum(n for i,n in read.cigartuples if i in {4,5})} 
                info={'strand':strand, 'chr':chrom,
                        'ID':read.query_name,'transcripts':{read.query_name:tr}}
                cov={}
                for read_id in tags['im'].split(','):
                    run=read_id[:read_id.find('/')]
                    runs.add(run)
                    cov[run] = cov.get(run, 0) + 1
                tr['coverage']=cov
                if 'SA' in tags: #part of a fusion gene   
                    if read.flag & 0x800: #secondary part
                        fusion.setdefault(read.query_name, []).append({'chr':chrom,'exons':tr['exons'],'cigar':tr['cigar']}) #cigar is required to sort out the order 
                    else: #primary part
                        fusion_primary[read.query_name]=tr
                if not read.flag & 0x800:
                    genes[chrom].add(Gene(exons[0][0],exons[-1][1],info))
        
        # make coverage a list
        runs=sorted(list(runs))
        for tree in genes.values():
            for g in tree:
                tr=next(iter(g.transcripts.values()))
                tr['coverage']=[tr['coverage'].get(run,0) for run in runs ]        
        #link secondary alignments (chimeric genes)
        for tid,primary in fusion_primary.items():
            primary['fusion']=fusion[tid]
        if n_secondary>0:
            log.info(f"skipped {n_secondary} secondary transcripts {n_secondary/total_reads*100}%")
        return genes,{'runs':runs}

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
            
        else:
            skipped.add(ls[2])
    if skipped:
        log.info('skipped the following categories: {}'.format(skipped))
    # sort the exons
    for tid in exons.keys():
        exons[tid].sort()
        
    missed_genes=[gid for gid in transcripts.keys() if gid not in gene_set]
    if missed_genes:        
        #log.debug('/n'.join(gid+str(tr) for gid, tr in missed_genes.items()))
        notfound=len(missed_genes)
        found=sum((len(t) for t in genes.values()) )
        log.warning('found specific gene information with category {} for {}/{} genes'.format(gene_categories,found, found+notfound))
    # add transcripts to genes
    for chrom in genes:
        for gene in genes[chrom]:
            g_id = gene.data['ID']
            tr=transcripts.get(g_id,  {g_id:{}})
            for t_id,tr_info in tr.items():
                try:
                    tr_info['exons']=exons[t_id]
                except KeyError:
                    # genes without transcripts get a single exons transcript
                    tr_info['exons']=[tuple(gene[:2])]
                gene.data.setdefault('transcripts', dict())[t_id]=tr_info
    return genes,dict()

def get_gff_chrom_dict(gff, chromosomes):
    # fetch chromosome ids
    chrom = {}
    for c in gff.contigs:
        #loggin.debug ("---"+c)
        for line in gff.fetch(c, 1, 2):
            if line[1] == "C":
                ls = line.split(sep="\t")
                if ls[2] == "region":
                    info = dict([pair.split("=")
                                    for pair in ls[8].split(";")])
                    if "chromosome" in info.keys():
                        if chromosomes is None or info["chromosome"] in chromosomes:
                            chrom[ls[0]] = info["chromosome"]
                        break
                        
        else:# no specific regions entrie - no aliases
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
            new_data={'chr':g.chrom, 'ID':f'{gene_prefix}{n_novel+offset}', 'strand':g.strand, 'transcripts':new_transcripts}
            genes.add(Gene(g.start,g.end,new_data ))
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

'''
def is_truncation(exons1, exons2):    
    relation=get_relation(exons1, exons2)
    if any(len(r) > 1 for r in relation):#
        log.debug('one exon from 1 has several corresponding exons from 2')
        return False  
    if not any(r for r in relation ): #
        log.debug('no relation')
        return False 
    first=next(iter([i for i,r in enumerate(relation) if r]))
    last=len(relation)-next(iter([i for i,r in enumerate(reversed(relation)) if r]))-1
    #if sum(1 for r in relation if r) != last-first+1:  # 
    if any(len(r)!=1 for r in relation[first:(last+1)]):
        log.debug('some exons in 1 do not have corresponding exons in 2')
        return False
    if first>0 and relation[first][0][0]>0:
        log.debug('start exons do not correspond to each other')
        return False # 
    if last < len(exons1)-1 and relation[last][0][0]<len(exons2)-1 :
        log.debug('last exons do not correspond to each other')
        return False # 
    if last!=first and ( #more than one exon
            (relation[first][0][1] & 1 ==0) or # 2nd splice junction must match
            (relation[last][0][1] & 2 == 0)): # fst splice junction must match
        log.debug('more than one exon, 2nd splice junction must match, fst splice junction must match')
        return False        
    if (relation[first][0][0]!=first and # one of the transcripts has extra exons at start
            (first==0 ) == (exons1[first][0] < exons2[relation[first][0][0]][0])): # check the begin of the shorter
        log.debug('first exon of trunkated version is actually longer')
        return False
    if (relation[last][0][0] != len(exons2)-1 and # one of the transcripts has extra exons at the end
            (last==len(exons1)-1) == (exons1[last][1] > exons2[relation[last][0][0]][1])):  #check the end of the shorter
        log.debug('last exon of trunkated version is actually longer')
        return False
    if last-first > 1 and any(relation[i][0][1]!=3 for i in range(first+1, last)): #
        log.debug('intermediate exons do not fit')
        return False
    log.debug('all filter passed')
    return relation, first, last
'''

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

def get_support(exons, ref_genes, chrom, is_reverse):        
    if chrom not in ref_genes.data:
        return None
    strand='-' if is_reverse else '+'
    ref_genes_ol = [g for g in ref_genes.data[chrom][exons[0][0]: exons[-1][1]] if g.strand==strand]
    #compute support for all transcripts of overlapping genes
    support_dict, novel_sj1, novel_sj2 = compute_support(ref_genes_ol, exons)
    #chose the best transcript
    try:
        # https://stackoverflow.com/questions/2474015/getting-the-index-of-the-returned-max-or-min-item-using-max-min-on-a-list
        best_idx = max(enumerate(zip(support_dict['sjIoU'], support_dict['exIoU'])), key=lambda x: x[1])[0]
    except ValueError:
        return None
    else:
        if support_dict['exIoU'][best_idx] <=0:
            return None
        support={'novel_sj1':novel_sj1, 'novel_sj2':novel_sj2}
        for n, v in support_dict.items():
            support[n]=v[best_idx]
        support["sType"]=get_alternative_splicing(support_dict, best_idx,exons, ref_genes_ol, is_reverse)
    return support

def compute_support(ref_genes, query_exons): #compute the support of all transcripts in ref_genes
    # get overlapping genes:
    query_len = sum([e[1]-e[0] for e in query_exons])
    query_nsj = len(query_exons)*2-2
    support = {k: list() for k in ['ref_gene_name','ref_gene_id', 'ref_transcript',
                                   'ref_tss', 'ref_pas', 'ref_len', 'ref_nSJ', 'exI', 'sjI']}    
    novel_sj1=[i+1 for i,e in enumerate(query_exons[1:]) if e[0] not in (ref_e[0] for g in ref_genes for tr in g.transcripts.values() for ref_e in tr['exons'])]
    novel_sj2=[i for i,e in enumerate(query_exons[:-1]) if e[1] not in (ref_e[1] for g in ref_genes for tr in g.transcripts.values() for ref_e in tr['exons'])]
    #novel_exons=[i for i,e in enumerate(query_exons[1:-1]) if e not in (ref_e) for g in ref_genes for tr in g.transcripts for ref_e in tr['exons'])]
    
    for gene in ref_genes:
        # tood: check strand??
        for trid,tr in gene.transcripts.items():
            db_exons=tr['exons']
            support['ref_gene_id'].append(gene.id)
            support['ref_gene_name'].append(gene.name)
            support['ref_transcript'].append(trid)
            support['ref_len'].append(sum([e[1]-e[0] for e in db_exons]))
            support['ref_nSJ'].append(len(db_exons)*2-2)
            support['ref_tss'].append(
                db_exons[0][0] if gene.data['strand'] == '+' else db_exons[-1][1])
            support['ref_pas'].append(
                db_exons[0][0] if gene.data['strand'] == '-' else db_exons[-1][1])
            sji, exi = get_intersects(db_exons, query_exons)
            support["sjI"].append(sji)
            support['exI'].append(exi)
    support['exIoU'] = [i/(query_len+db-i)
                        for i, db in zip(support['exI'],  support['ref_len'])]
    support['sjIoU'] = [i/(query_nsj+db-i) if (query_nsj+db-i) >
                        0 else 1 for i, db in zip(support['sjI'],  support['ref_nSJ'])]
    return support, novel_sj1, novel_sj2

def get_alternative_splicing(support_dict, best_idx,exons, ref_genes_ol, is_reverse):
    if support_dict["sjIoU"][best_idx] == 1:
        splice_type = ['splice_identical']
    elif support_dict["exIoU"][best_idx] == 0:
        splice_type = ['novel/unknown']
    else:
        try:
            g=next(iter(gene for gene in ref_genes_ol if gene.id==support_dict['ref_gene_id'][best_idx] ))
            ref_exons=next(iter(tr['exons'] for trid,tr in g.transcripts.items() if trid == support_dict['ref_transcript'][best_idx]))
        except StopIteration:
            log.error('cannot find the transcript {}:{}-- this should never happen'.format(*[support_dict[v][best_idx] for v in ['ref_gene_name', 'ref_transcript']]))
            raise ValueError
        splice_type = get_splice_type( ref_exons, exons, is_reverse)
        fusion=get_fusion(support_dict, ref_genes_ol)
        if len(fusion)>0:
            splice_type['fusion_gene']=fusion
    return splice_type

def get_fusion(support_dict, ref_genes_ol, min_iou=.5):
        covered_gene_ids = {id for id, iou in zip(
            support_dict['ref_gene_id'], support_dict['sjIoU']) if iou >= min_iou}
        covered_genes = {
            g for g in ref_genes_ol if g.id in covered_gene_ids}        
        fusion=[]
        if len(covered_genes) > 1:
            # remove overlapping genes
            for g1, g2 in combinations(covered_genes, 2):
                if not overlap(g1[:2], g2[:2]):
                    fusion.append('{}|{}'.format(g1.name, g2.name))
        return fusion
    

def get_splice_type(ref, alt, is_reversed=False):
    if len(ref) == 0:
        return(['novel/unknown'])
    if len(alt) ==1:
        return['unspliced_fragment']
    types = ['alternative_donor', 'alternative_acceptor', 'alternative_promoter', 'alternative_polyA',
             'truncated5', 'truncated3', 'exon_skipping', 'novel_exon', 'gapped_exon', 'retained_intron']
    types = {t: [] for t in types}
    relation = get_relation(ref, alt)
    # alt exons that do not overlap any ref
    first = next(i for i, rel in enumerate(relation) if rel) #first ref exon that overlaps an alt exon
    last=len(relation)-next(iter([i for i,r in enumerate(reversed(relation)) if r]))-1
    novel = set(range(len(alt)))-{r[0] for x in relation for r in x}
    if 0 in novel: #frist exon completely new (wrt troi)
        types['alternative_polyA' if is_reversed else 'alternative_promoter'].append(
            '{}-{}'.format(*alt[0]))
    if len(alt)-1 in novel:
        types['alternative_promoter' if is_reversed else 'alternative_polyA'].append(
            '{}-{}'.format(*alt[-1]))
    if len(novel - {0, len(alt)-1}) > 0:
        # todo: distinguish completely novel from novel for this ref_transcript
        types['novel_exon'] += ['{}-{}'.format(*e)
                                for i, e in enumerate(alt) if i in novel-{0, len(alt)-1}]
    # if 0 not in novel and not types['novel_exon'] and len(relation[0]) == 0:
    all_match = (len(novel) == 0 and #no novel alt exons
        all([len(rel) == 1 for rel in relation[first:(last+1)]]) and #each alt exon corresponds to exactly one ref exon
        all([rel[0][1] == 3 for rel in relation[(first+1):last]]) and #all but the first and last alt exons are the same
        (relation[first][0][1] & 1) and #the first exon has matching second splice site
        (relation[last][0][1] & 2)) #the last exon has matching first splice site
    if first > 0 and all_match :
        if alt[0][0] >= ref[first][0]: #todo: check that alt[0][0]> ref[first][0]
            types['truncated3' if is_reversed else 'truncated5'].append(
                '{}-{}'.format(*alt[0]))
        else:
            types['alternative_polyA' if is_reversed else 'alternative_promoter'].append(
                '{}-{}'.format(*alt[0]))
    # if len(alt)-1 not in novel and not types['novel_exon'] and len(relation[-1]) == 0:
    if last < len(relation)-1 and all_match:
        if alt[-1][1] <= ref[last][1]:
            types['truncated5' if is_reversed else 'truncated3'].append(
                '{}-{}'.format(*alt[-1]))
        else:
            types['alternative_promoter' if is_reversed else 'alternative_polyA'].append(
                '{}-{}'.format(*alt[-1]))
    for i in range(first, last+1):
        rel = relation[i]
        if len(rel) > 1:  # more than one alt exon for a ref exon
            types['gapped_exon'].append(
                 "~".join(['{}-{}'.format(*alt[alt_i]) for alt_i, splice in rel]))
        if rel and i > first and relation[i-1] and rel[0][0] == relation[i-1][-1][0]:
            types['retained_intron'].append('{}-{}'.format(*alt[rel[0][0]]))
        else:
            # fst splice site does not correspond to any alt exon
            if rel and rel[0][1] & 2 == 0 and i > first:
                delta = alt[rel[0][0]][0]-ref[i][0]
                types['alternative_acceptor' if is_reversed else 'alternative_donor'].append(
                    (alt[rel[0][0]][0], delta))
        if rel and i < last and rel[-1][1] & 1 == 0 and i < last and not(relation[i+1] and rel[-1][0] == relation[i+1][0][0]):
                delta = alt[rel[-1][0]][1]-ref[i][1]
                types['alternative_donor' if is_reversed else 'alternative_acceptor'].append(
                    (alt[rel[-1][0]][1], delta))
        if not rel and i > first and i < last:  # exon skipping
            pre=next(((j,relation[j][-1]) for j in reversed(range(i)) if relation[j]),[]) #first ref exon that overlaps an alt exon
            suc=next(((j,relation[j][0]) for j in range(i+1, len(relation)) if relation[j]),[]) #first ref exon that overlaps an alt exon
            # only predecessing as well as successing splice site is identical
            #if pre and suc and pre[1][1] & 1 and suc[1][1] & 2:
            types['exon_skipping'].append(
                    '{}~{}'.format(alt[pre[1][0]][1], alt[suc[1][0]][0]))
    # return only types that are present
    return {k: v for k, v in types.items() if v}
