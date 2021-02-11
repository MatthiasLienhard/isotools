import os
import pickle
import logging
from ._transcriptome_io import import_gtf_transcripts, import_gff_transcripts
from .gene import Gene
from intervaltree import IntervalTree #, Interval
import pandas as pd
logger=logging.getLogger('isotools')

# as this class has diverse functionality, its split among:
# transcriptome.py (this file- initialization and user level basic functions)
# _transcriptome_io.py (input/output primary data files/tables)
# _transcriptome_stats.py (statistical methods)
# _trnascriptome_plots.py (plots)
# _transcriptome_filter.py (gene/transcript iteration and filtering)

class Transcriptome:
    '''Container class for genes
    '' contains a dict of interval trees for the genes, each containing the splice graph'''
    #####initialization and save/restore data
    def __new__(cls, pickle_file=None,**kwargs):
        if pickle_file is not None:
            obj=cls.load(pickle_file)
        else:
            obj=super().__new__(cls)
        return obj

    def __init__(self, pickle_file=None,**kwargs ):     
        if 'data' in kwargs:
            self.data,self.infos, self.chimeric=kwargs['data'],kwargs.get('infos',dict()),kwargs.get('chimeric',{})
            assert 'reference_file' in self.infos 
            self.make_index()
        
            
    
    @classmethod
    def from_reference(cls, reference_file, file_format='auto',**kwargs):
        tr=cls.__new__(cls)
        tr.infos={'reference_file':reference_file}
        tr.chimeric={}
        if file_format=='auto':        
            file_format=os.path.splitext(reference_file)[1].lstrip('.')
            if file_format=='gz':
                file_format=os.path.splitext(reference_file[:-3])[1].lstrip('.')
        logger.info(f'importing reference from {file_format} file {reference_file}')
        if file_format == 'gtf':
            tr.data= import_gtf_transcripts(reference_file,tr,  **kwargs)
        elif file_format in ('gff', 'gff3'):            
            tr.data= import_gff_transcripts(reference_file,tr,  **kwargs)
        elif file_format == 'pkl':
            tr= pickle.load(open(reference_file, 'rb'))            
            if [k for k in tr.infos if k!='reference_file']:
                logger.warning('the pickle file seems to contain additional expression information... extracting refrence')
                tr=tr._extract_reference()
        else: logger.error('unknown file format %s',file_format)
        return tr

    @classmethod
    def load(cls, pickle_file):
        'restores the information of a transcriptome from a pickle file'
        logger.info('loading transcriptome from '+pickle_file)
        return pickle.load(open(pickle_file, 'rb'))
        

    def save_reference(self, fn=None):    
        'saves the reference information of a transcriptome in a pickle file'
        if fn is None:
            fn=self.infos['reference_file']+'.isotools.pkl'
        logger.info('saving reference to '+fn)       
        ref_tr=self._extract_reference() 
        pickle.dump(ref_tr, open(fn, 'wb'))

    def _extract_reference(self):
        if not [k for k in self.infos if k!='reference_file']:
            return self #only reference info - assume that self.data only contains reference data
        #make a new transcriptome
        ref_tr=type(self)(data={} ,infos={'reference_file':self.infos['reference_file']})
        # extract the reference genes and link them to the new ref_tr
        for chrom,tree in self.data.items():
            ref_tr.data[chrom]=IntervalTree(Gene(g.start,g.end,{k:g.data[k] for k in Gene.required_infos+['reference']}, ref_tr) for g in tree if g.is_annotated)
        return ref_tr

    def save(self, fn=None):
        'saves the information of a transcriptome (including reference) in a pickle file'
        if fn is None:
            fn=self.infos['out_file_name']+'.isotools.pkl' #key error if not set
        logger.info('saving transcriptome to '+fn)
        pickle.dump(self, open(fn, 'wb'))
    
    def make_index(self):
        'updates the index used for __getitem__, e.g. the [] operator'
        idx=dict()
        for g in self:
            if g.id in idx: # at least id should be unique - maybe raise exception?
                logger.warn(f'{g.id} seems to be ambigous: {str(self[g.id])} vs {str(g)}')
            idx[g.name] = g
            idx[g.id]=g
        self._idx=idx
 
    ##### basic user level functionality
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

    def _get_sample_idx(self, group_column='name'):
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

    ### IO: load new data from primary data files
    from ._transcriptome_io import add_sample_from_bam,remove_samples,add_short_read_coverage, collapse_immune_genes

    ### IO: utility functions
    from ._transcriptome_io import _add_sample_transcript, _add_novel_genes, _get_intersects, _add_chimeric
    
    ### IO: output data as tables or other human readable format
    from ._transcriptome_io import gene_table, transcript_table,chimeric_table,write_gtf

    ### filtering functionality and iterators
    from ._transcriptome_filter import add_biases, add_filter,iter_genes,iter_transcripts,iter_ref_transcripts

    ### statistic: differential splicing, embedding
    from ._transcriptome_stats import altsplice_test,splice_dependence_test, embedding

    # statistic: summary tables (can be used as input to plot_bar / plot_dist)
    from ._transcriptome_stats import altsplice_stats,filter_stats,transcript_length_hist,transcript_coverage_hist,transcripts_per_gene_hist,exons_per_transcript_hist,downstream_a_hist, direct_repeat_hist
