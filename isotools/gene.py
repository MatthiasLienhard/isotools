from intervaltree import Interval
from Bio.Seq import Seq
from itertools import combinations
import numpy as np
import copy
from .splice_graph import SpliceGraph
from ._utils import splice_identical
import logging

class Gene(Interval):
    'this class stores all gene information and transcripts. It is derived from intervaltree.Interval'
    required_infos=['ID','name', 'chr', 'strand']

    def __new__(cls,begin, end, data, transcriptome):
        return super().__new__(cls,begin, end, data) #required as Interval (and Gene) is immutable

    def __init__(self,begin, end, data, transcriptome):
        self._transcriptome=transcriptome

    def __str__(self):
        return('Gene {} {}({}), {} reference transcripts, {} expressed transcripts'.format(self.name, self.region, self.strand, self.n_ref_transcripts, self.n_transcripts))
    
    def __repr__(self):
        return object.__repr__(self)

    def add_filter(self, gene_filter,transcript_filter,ref_transcript_filter):   
        'add the filter flags'
        self.data['filter']=[name for name,fun in gene_filter.items() if  fun(**self.data)]
        for tr in self.transcripts:
            tr['filter']=[name for name,fun in transcript_filter.items() if fun(**tr)]
        for tr in self.ref_transcripts:
            tr['filter']=[name for name,fun in ref_transcript_filter.items() if fun(**tr)]

    def filter_transcripts(self, include, remove):
        ''' iterator over the transcripts of the gene. Filtering implemented by passing lists of flags to the parameters:
        include: transcripts must have at least one of the flags
        remove: transcripts must not have one of the flags'''
        for i,tr in enumerate(self.transcripts):
            if not include or any(f in tr['filter'] for f in include):
                if not remove or all(f not in tr['filter'] for f in remove):
                    yield i,tr 

    def filter_ref_transcripts(self, include, remove):
        ''' iterator over the refernce transcripts of the gene. Filtering implemented by passing lists of flags to the parameters:
        include: transcripts must have at least one of the flags
        remove: transcripts must not have one of the flags'''
        for i,tr in enumerate(self.ref_transcripts):
            if not include or any(f in tr['filter'] for f in include):
                if not remove or all(f not in tr['filter'] for f in remove):
                    yield i,tr
        
    def correct_fuzzy_junctions(self,tr, size, modify=True):
        'this corrects for splicing shifts (e.g. same difference compared to reference annotaiton at both donor and acceptor) presumably caused by ambigous alignments'
        exons=tr['exons']
        shifts=self.ref_splice_graph.fuzzy_junction(exons,size)
        if shifts and modify:
            for i,sh in shifts.items():
                exons[i][1]+=sh
                exons[i+1][0]+=sh
                
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
        #todo: this is highly redundant - should be done on splice graph
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
        a_content={}
        for tr in (t for tL in (self.transcripts, self.ref_transcripts) for t in tL):
            if self.strand=='+':
                pos=tr['exons'][-1][1]
            else:
                pos=tr['exons'][0][0]-length    
            if pos not in a_content:                
                seq=genome_fh.fetch(self.chrom, max(0,pos), pos+length) 
                if self.strand=='+':
                    a_content[pos]=seq.upper().count('A')/length
                else:
                    a_content[pos]=seq.upper().count('T')/length
            tr['downstream_A_content']=a_content[pos]

    def add_truncations(self): 
        'check for truncations and save indices of truncated transcripts'
        for trid, containers in self.splice_graph.find_truncations().items():
            self.transcripts[trid]['truncated']=containers # list of transcript ids containing trid
        #todo: add this information how much is truncated: transcripts[trid1]['truncated']=(trid2,delta5, delta3)
        
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
            logging.error(self.data)
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
            return []

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
            self.data['reference']['splice_graph']=SpliceGraph(exons)
        return self.data['reference']['splice_graph']
        
    @property
    def splice_graph(self):
        if 'splice_graph' not in self.data or self.data['splice_graph'] is None:        
            samples=self._transcriptome.samples            
            weights=np.zeros((len(samples), self.n_transcripts))
            for i,sa in enumerate(samples):
                for j,tr in enumerate(self.transcripts):
                    try:
                        weights[i,j]=sum(cov for cov in tr['samples'][sa]['range'].values())
                    except KeyError:
                        pass
            exons=[tr['exons'] for tr in self.transcripts]
            self.data['splice_graph']=SpliceGraph(exons, samples, weights)
        return self.data['splice_graph']

    def __copy__(self):
        return Gene(self.start, self.end, self.data, self._transcriptome)        
        

    def __deepcopy__(self, memo):
        return Gene(self.start, self.end, copy.deepcopy(self.data, memo), self._transcriptome)
        

    def __reduce__(self):
        return Gene, (self.start, self.end, self.data, self._transcriptome)  

    def copy(self):
        return self.__copy__()


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

