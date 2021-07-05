from intervaltree import Interval
from Bio.Seq import Seq
import numpy as np
import copy
from .splice_graph import SegmentGraph
from .short_read import Coverage

import logging
logger=logging.getLogger('isotools')

def _eval_filter_fun(fun,name,args):
    '''Decorator for the filter functions, which are lambdas and thus cannot have normal decorators.
    On exceptions the provided parameters are reported. This is helpfull for debugging.'''
    try:
        return fun(**args)
    except Exception as e:
        logger.error('error when evaluating filter %s with arguments %s: %s',name,str(args), str(e))
        raise           #either stop evaluation 
        #return False   #or continue


class Gene(Interval):
    'This class stores all gene information and transcripts. It is derived from intervaltree.Interval.'
    required_infos=['ID','name', 'chr', 'strand']

    ###initialization
    def __new__(cls,begin, end, data, transcriptome):
        return super().__new__(cls,begin, end, data) #required as Interval (and Gene) is immutable

    def __init__(self,begin, end, data, transcriptome):
        self._transcriptome=transcriptome

    def __str__(self):
        return('Gene {} {}({}), {} reference transcripts, {} expressed transcripts'.format(self.name, self.region, self.strand, self.n_ref_transcripts, self.n_transcripts))
    
    def __repr__(self):
        return object.__repr__(self)

    from ._gene_plots import sashimi_plot,gene_track,sashimi_plot_short_reads, sashimi_figure

    def short_reads(self,idx):
        '''Returns the short read coverage profile for a short read sample. 
        
        :param idx: The index of the short read sample. '''

        try:
            return self.data['short_reads'][idx]
        except (KeyError, IndexError):
            srdf=self._transcriptome.infos['short_reads'] #raises key_error if no short reads added
            self.data.setdefault('short_reads',[])
            for i in range(len(self.data['short_reads']),len(srdf)):
                self.data['short_reads'].append(Coverage.from_bam(srdf.file[i],self))
        return self.data['short_reads'][idx]

    def add_filter(self, gene_filter,transcript_filter,ref_transcript_filter):   
        '''Adds the filter flags to the gene, transcripts and reference transcripts.
        
        This evaluates the provided expressions on the gene, transcripts and reference transcripts and adds the appropriate filter flags.

        :param gene_filter: Flag expressions for the gene level information. 
        :param transcript_filter: Flag expressions to be evalutated on the transcripts. 
        :param ref_transcript_filter: Flag expressions to be evalutated on the reference transcripts. '''
        self.data['filter']=[name for name,fun in gene_filter.items() if  _eval_filter_fun(fun,name,self.data)]
        for tr in self.transcripts:
            tr['filter']=[name for name,fun in transcript_filter.items() if _eval_filter_fun(fun,name,tr)]
        for tr in self.ref_transcripts:
            tr['filter']=[name for name,fun in ref_transcript_filter.items() if _eval_filter_fun(fun,name,tr)]

    def filter_transcripts(self, include, remove, anno_include, anno_remove, mincoverage=None, maxcoverage=None):
        ''' Iterator over the transcripts of the gene. 
        
        Transcrips are specified by lists of flags submitted to the parameters.
        
        :param include: transcripts must have at least one of the flags
        :param remove: transcripts must not have one of the flags'''
        for i,tr in enumerate(self.transcripts):
            if mincoverage and self.coverage[:,i].sum()< mincoverage:
                continue
            if maxcoverage and self.coverage[:,i].sum()> maxcoverage:
                continue
            if not include or any(f in tr['filter'] for f in include):
                if not anno_include or any(f in tr['annotation'][1] for f in anno_include):
                    if not remove or all(f not in tr['filter'] for f in remove):
                        if not anno_remove or all(f not in tr['annotation'][1] for f in anno_remove):
                            yield i,tr 

    def filter_ref_transcripts(self, include, remove):
        ''' Iterator over the refernce transcripts of the gene. 
        
        Transcrips are specified by lists of flags submitted to the parameters.

        :param include: Transcripts must have at least one of the flags.
        :param remove: Transcripts must not have one of the flags.'''
        for i,tr in enumerate(self.ref_transcripts):
            if not include or any(f in tr['filter'] for f in include):
                if not remove or all(f not in tr['filter'] for f in remove):
                    yield i,tr
        
    def correct_fuzzy_junctions(self,trid, size, modify=True):
        '''Corrects for splicing shifts.
        
         This function looks for "shifted junctions", e.g. same difference compared to reference annotaiton at both donor and acceptor) 
         presumably caused by ambigous alignments. In these cases the positions are adapted to the reference position (if modify is set). 
         
         :param trid: The index of the transcript to be checked.
         :param size: The maximum shift to be corrected. 
         :param modify: If set, the exon positions are corrected according to the reference.'''

        exons=trid['exons']
        shifts=self.ref_segment_graph.fuzzy_junction(exons,size)
        if shifts and modify:
            for i,sh in shifts.items():
                exons[i][1]+=sh
                exons[i+1][0]+=sh
                
        return shifts
 
    def _to_gtf(self,trids, source='isoseq', use_gene_name=False):
        '''Creates the gtf lines of the gene as strings.'''

        info={'gene_id':self.name if use_gene_name else self.id}
        lines=[]
        starts=[]
        ends=[]
        for i in trids:
            tr=self.transcripts[i]
            info['transcript_id']='{}.{}'.format(info['gene_id'],i)
            #todo: print only relevant infos
            starts.append(tr['exons'][0][0]+1)
            ends.append(tr['exons'][-1][1])
            lines.append((self.chrom, source, 'transcript', tr['exons'][0][0]+1, tr['exons'][-1][1], '.',self.strand, '.', '; '.join(f'{k} "{v}"' for k,v in info.items() if k != 'exon_id')))
            for enr, pos in enumerate(tr['exons']):
                info['exon_id']='{}.{}_{}'.format(self.id,i, enr)
                lines.append((self.chrom, source, 'exon', pos[0]+1, pos[1], '.',self.strand, '.', '; '.join(f'{k} "{v}"' for k,v in info.items())))
        if lines:
            #add gene line
            lines.append((self.chrom, source, 'gene', min(starts), max(ends), '.',self.strand, '.', '{} "{}"'.format('gene_id',info['gene_id'])))
        return lines
 
    def add_noncanonical_splicing(self, genome_fh):
        '''Add information on noncanonical splicing. 
        
        For all transcripts of the gene, scan for noncanonical (i.e. not GT-AG) splice sites.
        If noncanonical splice sites are present, the corresponding intron index (in genomic orientation) and the sequence
        i.e. the dinucleotides of donor and aceptor as XX-YY string are stored in the "noncannoncical_splicing" field of the transcript dicts.
        True noncanonical splicing is rare, thus it might indicate technical artifacts (template switching, missalignment, ...)
        
        :param genome_fh: A file handle of the genome fasta file.'''
        ss_seq={}
        for tr in self.transcripts:
            pos=[(tr['exons'][i][1], tr['exons'][i+1][0]-2) for i in range(len(tr['exons'])-1)]
            new_ss_seq={   site:genome_fh.fetch(self.chrom, site , site+2).upper() for intron in pos for site in intron if site not in ss_seq  }
            if new_ss_seq:
                ss_seq.update(new_ss_seq)

            if self.strand=='+':
                sj_seq=[ss_seq[d]+ss_seq[a] for d,a in pos]
            else:
                sj_seq=[str(Seq(ss_seq[d]+ss_seq[a]).reverse_complement()) for d,a in pos]
            
            nc=[(i,seq) for i,seq in enumerate(sj_seq) if seq != 'GTAG']
            if nc:
                tr['noncanonical_splicing']=nc

    def add_direct_repeat_len(self, genome_fh,delta=10):
        '''Computes direct repeat length. 
        
        This function counts the number of consequtive equal bases at donor and acceptor sites of the splice junctions. 
        This information is stored in the "direct_repeat_len" filed of the transcript dictionaries. 
        Direct repeats longer than expected by chance indicate template switching.
        
        :param genome_fh: The file handle to the genome fasta.
        :param delta: The maximum length of direct repeats that can be found. '''
        
        intron_seq={}
        score={}
        for tr in self.transcripts:
            for intron in ((tr['exons'][i][1], tr['exons'][i+1][0]) for i in range(len(tr['exons'])-1)):
                for pos in intron:
                    intron_seq.setdefault(pos, genome_fh.fetch(self.chrom, pos-delta, pos+delta))
                if intron not in score:
                    #align=[pairwise2.align.globalms(seq5, seq3, 1,-1, -1, 0) for seq5, seq3 in intron_seq] 
                    #align=[pairwise2.align.globalxs(seq5, seq3,-1,-1, score_only=True) for seq5, seq3 in intron_seq] # number of equal bases at splice site
                    align=[a==b for a,b in zip(intron_seq[intron[0]],intron_seq[intron[1]])]
                    score[intron]=0
                    for a in align[delta:]: #find the runlength at the middle (delta)
                        if not a:
                            break
                        score[intron]+=1
                    for a in reversed(align[:delta]):
                        if not a:
                            break
                        score[intron]+=1
        for tr in self.transcripts:
            tr['direct_repeat_len']=[min(score[intron], delta) for intron in ((tr['exons'][i][1], tr['exons'][i+1][0]) for i in range(len(tr['exons'])-1))]
    
    def add_threeprime_a_content(self, genome_fh, length=30):
        '''Adds the information of the genomic A content downstream the transcript. 

        High values of genomic A content indicate internal priming and hence genomic origin of the LRTS read. 
        This function populates the 'downstream_A_content' field of the transcript dictionaries.

        :param geneome_fh: A file handle for the indexed genome fasta file. 
        :param length: The length of the downstream region to be considered. 
        '''
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

    def add_fragments(self): 
        '''Checks for transcripts that are fully contained in other transcripts.
        
        Transcripts that are fully contained in other transcripts are potential truncations. 
        This function populates the 'fragment' filed of the transcript dictionaries with the indices of the containing transcripts, 
        and the exon ids that match the first and last exons.'''
            
        for trid, containers in self.segment_graph.find_fragments().items():
            self.transcripts[trid]['fragments']=containers # list of (containing transcript id, first 5' exons, first 3'exons)
       
        
    def coding_len(self, trid):
        '''Returns length of 5\'UTR, coding sequence and 3\'UTR.
        
        :param trid: The transcript index for which the coding length is requested. '''

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

    def get_infos(self,trid, keys):
        '''Returns the transcript information specified in "keys" as a list.'''
        return [value for k in keys for value in self._get_info(trid,k)]


    def _get_info(self ,trid, key):
        #returns tuples (as some keys return multiple values)
        if key=='length':
            return sum((e-b for b,e in self.transcripts[trid]['exons'])),
        if key=='filter':
            if self.transcripts[trid]['filter']:
                return ','.join(self.transcripts[trid]['filter']),
            else:
                return 'PASS',
        elif key=='n_exons':
            return len(self.transcripts[trid]['exons']),
        elif key=='exon_starts':
            return ','.join(str(e[0]) for e in self.transcripts[trid]['exons']),
        elif key=='exon_ends':
            return ','.join(str(e[1]) for e in self.transcripts[trid]['exons']),
        elif key=='annotation':
            #sel=['sj_i','base_i', 'as']
            if 'annotation' not in self.transcripts[trid]:
                return ('NA',)*2
            nov_class, subcat=self.transcripts[trid]['annotation']
            subcat_string=';'.join(k if v is None else '{}:{}'.format(k,v) for k,v in subcat.items())
            return nov_class,subcat_string
        elif key=='coverage':
            return self.coverage[:,trid]
        elif key in self.transcripts[trid]:
            val=self.transcripts[trid][key]
            try:
                iter(val)
            except TypeError:
                return val, #atomic (e.g. numeric)
            else:
                return str(val), #iterables get converted to string
        return 'NA',

    def _set_coverage(self, force=False):
        samples=self._transcriptome.samples
        cov=np.zeros((len(samples), self.n_transcripts), dtype=int)
        if not force: #keep the segment graph if no new transcripts
            known=self.data.get('coverage', None)            
            if known is not None and known.shape[1]==self.n_transcripts:
                if known.shape==cov.shape: 
                    return
                cov[:known.shape[0],:]=known
                for i in range(known.shape[0],len(samples)):
                    for j,tr in enumerate(self.transcripts):
                        cov[i,j]=tr['coverage'].get(samples[i],0)
                self.data['coverage']=cov
                return
        for i,sa in enumerate(samples):
            for j,tr in enumerate(self.transcripts):
                    cov[i,j]=tr['coverage'].get(sa,0)
        self.data['coverage']=cov
        self.data['segment_graph']=None


    @property
    def coverage(self):
        '''Returns the transcript coverage.

        Coverage is returned as a numpy array, with samples in columns and transcript isoforms in the rows.'''
        cov=self.data.get('coverage', None)
        if cov is not None:
            return cov
        self._set_coverage()
        return self.data['coverage']

    @property
    def gene_coverage(self):
        '''Returns the gene coverage.
        
        Total Coverage of the gene for each sample.'''

        return self.coverage.sum(1)

    @property
    def chrom(self):
        '''Returns the genes chromosome.'''
        return self.data['chr']
        
    @property
    def start(self): #alias for begin
        return self.begin
    
    @property
    def region(self):
        '''Returns the region of the gene as a string in the format "chr:start-end".'''
        try:
            return '{}:{}-{}'.format(self.chrom,self.start, self.end )
        except KeyError:
            raise

    @property
    def id(self):
        '''Returns the gene id'''
        try:
            return self.data['ID']    
        except KeyError:
            logger.error(self.data)
            raise
            
    @property
    def name(self):
        '''Returns the gene name'''
        try:
            return self.data['name']
        except KeyError:
            return self.id #e.g. novel genes do not have a name (but id)
    
    @property
    def is_annotated(self):
        '''Returns "True" iff reference annotation is present for the gene.'''
        return 'reference' in self.data
    
    @property
    def is_expressed(self):
        '''Returns "True" iff gene is covered by at least one long read in at least one sample.'''
        return bool(self.transcripts)

    @property
    def strand(self):
        '''Returns the strand of the gene, e.g. "+" or "-"'''
        return self.data['strand']

    @property
    def transcripts(self):
        '''Returns the list of transcripts of the gene, as found by LRTS.'''
        try:
            return self.data['transcripts']
        except KeyError:
            return []

    @property
    def ref_transcripts(self):
        '''Returns the list of reference transcripts of the gene.'''
        try:
            return self.data['reference']['transcripts']
        except KeyError:
            return []

    @property
    def n_transcripts(self):
        '''Returns number of transcripts of the gene, as found by LRTS.'''
        return len(self.transcripts)
    
    @property
    def n_ref_transcripts(self):
        '''Returns number of reference transcripts of the gene.'''
        return len(self.ref_transcripts)

    @property
    def ref_segment_graph(self): #raises key error if not self.is_annotated
        '''Returns the segment graph of the reference transcripts for the gene'''

        assert self.is_annotated, "reference segment graph requested on novel gene"
        if 'segment_graph' not in self.data['reference'] or self.data['reference']['segment_graph'] is None:
            exons=[tr['exons'] for tr in self.ref_transcripts]
            self.data['reference']['segment_graph']=SegmentGraph(exons, self.strand)
        return self.data['reference']['segment_graph']
        
    @property
    def segment_graph(self):
        '''Returns the segment graph of the LRTS transcripts for the gene'''
        if 'segment_graph' not in self.data or self.data['segment_graph'] is None:        
            exons=[tr['exons'] for tr in self.transcripts]
            self.data['segment_graph']=SegmentGraph(exons, self.strand)
        return self.data['segment_graph']

    def __copy__(self):
        return Gene(self.start, self.end, self.data, self._transcriptome)        

    def __deepcopy__(self, memo): #does not copy _transcriptome!
        return Gene(self.start, self.end, copy.deepcopy(self.data, memo), self._transcriptome)

    def __reduce__(self):
        return Gene, (self.start, self.end, self.data, self._transcriptome)  

    def copy(self):
        'Returns a shallow copy of self.'
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

