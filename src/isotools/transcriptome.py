import os
import pickle
import logging
import re
from intervaltree import IntervalTree  # , Interval
import pandas as pd
from ._transcriptome_io import import_ref_transcripts
from .gene import Gene
from ._transcriptome_filter import DEFAULT_GENE_FILTER, DEFAULT_TRANSCRIPT_FILTER, DEFAULT_REF_TRANSCRIPT_FILTER, ANNOTATION_VOCABULARY, SPLICE_CATEGORY
from . import __version__
logger = logging.getLogger('isotools')


# as this class has diverse functionality, its split among:
# transcriptome.py (this file- initialization and user level basic functions)
# _transcriptome_io.py (input/output primary data files/tables)
# _transcriptome_stats.py (statistical methods)
# _trnascriptome_plots.py (plots)
# _transcriptome_filter.py (gene/transcript iteration and filtering)

class Transcriptome:
    '''Contains sequencing data and annotation for Long Read Transcriptome Sequencing (LRTS) Experiments.
    '''
    # initialization and save/restore data

    def __init__(self, **kwargs):
        '''Constructor method'''
        if 'data' in kwargs:
            self.data = kwargs['data']
            self.infos = kwargs.get('infos', {})
            self.chimeric = kwargs.get('chimeric', {})
            self.filter = kwargs.get('filter', {})
            assert 'reference_file' in self.infos
            self.make_index()

    @classmethod
    def from_reference(cls, reference_file, file_format='auto', **kwargs):
        '''Creates a Transcriptome object by importing reference annotation.

        :param reference_file: Reference file in gff3 format or pickle file to restore previously imported annotation
        :type reference_file: str
        :param file_format: Specify the file format of the provided reference_file.
            If set to "auto" the file type is infrered from the extension.
        :param chromosome: If reference file is gtf/gff, restrict import on specified chromosomes
        :param infer_transcrpts: If reference file is gtf, genes and transcripts are infered from "exon" entries, no specific transcript '''

        if file_format == 'auto':
            file_format = os.path.splitext(reference_file)[1].lstrip('.')
            if file_format == 'gz':
                file_format = os.path.splitext(reference_file[:-3])[1].lstrip('.')
        if file_format in ('gff', 'gff3', 'gtf'):
            logger.info('importing reference from %s file %s', file_format, reference_file)
            tr = cls()
            tr.chimeric = {}
            tr.data = import_ref_transcripts(reference_file, tr, file_format,  **kwargs)
            tr.infos = {'reference_file': reference_file, 'isotools_version': __version__}
            tr.filter = {'gene': DEFAULT_GENE_FILTER.copy(),
                         'transcript': DEFAULT_TRANSCRIPT_FILTER.copy(),
                         'reference': DEFAULT_REF_TRANSCRIPT_FILTER.copy()}
            for subcat in ANNOTATION_VOCABULARY:
                tag = '_'.join(re.findall(r'\b\w+\b', subcat)).upper()
                if tag[0].isdigit():
                    tag = '_'+tag
                tr.filter['transcript'][tag] = f'"{subcat}" in annotation[1]'
            for i, cat in enumerate(SPLICE_CATEGORY):
                tr.filter['transcript'][cat] = f'annotation[0]=={i}'

        elif file_format == 'pkl':
            # warn if kwargs are specified: kwargs are ignored
            if kwargs:
                logger.warning("The following parameters are ignored when loading reference from pkl: %s", ", ".join(kwargs))
            tr = cls.load(reference_file)
            if 'sample_table' in tr.infos:
                logger.warning('the pickle file seems to contain sample information... extracting refrence')
                tr = tr._extract_reference()
        else:
            logger.error('unknown file format %s', file_format)
        return tr

    def save(self, pickle_file=None):
        '''Saves transcriptome information (including reference) in a pickle file.

        :param pickle_file: Filename to save data'''
        if pickle_file is None:
            pickle_file = self.infos['out_file_name']+'.isotools.pkl'  # key error if not set
        logger.info('saving transcriptome to %s', pickle_file)
        pickle.dump(self, open(pickle_file, 'wb'))

    @classmethod
    def load(cls, pickle_file):
        '''Restores transcriptome information from a pickle file.

        :param pickle_file: Filename to restore data'''

        logger.info('loading transcriptome from %s', pickle_file)
        tr = pickle.load(open(pickle_file, 'rb'))
        pickled_version = tr.infos.get('isotools_version', '<0.2.6')
        if pickled_version != __version__:
            logger.warning('This is isotools version %s, but data has been pickled with version %s, which may be incompatible', __version__, pickled_version)
        return tr

    def save_reference(self, pickle_file=None):
        '''Saves the reference information of a transcriptome in a pickle file.

        :param pickle_file: Filename to save data'''
        if pickle_file is None:
            pickle_file = self.infos['reference_file']+'.isotools.pkl'
        logger.info('saving reference to %s', pickle_file)
        ref_tr = self._extract_reference()
        pickle.dump(ref_tr, open(pickle_file, 'wb'))

    def _extract_reference(self):
        if not 'sample_table' not in self.infos:
            return self  # only reference info - assume that self.data only contains reference data
        # make a new transcriptome
        ref_info = {k: v for k, v in self.infos.items() if k in ['reference_file', 'isotools_version']}
        ref_tr = type(self)(data={}, infos=ref_info, filter=self.filter)
        # extract the reference genes and link them to the new ref_tr
        keep = {'ID', 'chr', 'strand', 'name', 'reference'}  # no coverage, segment_graph, transcripts
        for chrom, tree in self.data.items():
            ref_tr.data[chrom] = IntervalTree(Gene(g.start, g.end, {k: v
                                              for k, v in g.data.items() if k in keep}, ref_tr) for g in tree if g.is_annotated)
        return ref_tr

    def make_index(self):
        '''Updates the index of gene names and ids (e.g. used by the the [] operator).'''
        idx = dict()
        for g in self:
            if g.id in idx:  # at least id should be unique - maybe raise exception?
                logger.warning('%s seems to be ambigous: %s vs %s', g.id, str(idx[g.id]), str(g))
            idx[g.name] = g
            idx[g.id] = g
        self._idx = idx

    # basic user level functionality
    def __getitem__(self, key):
        '''
        Syntax: self[key]

        :param key: May either be the gene name or the gene id
        :return: The gene specified by key.'''
        return self._idx[key]

    def __len__(self):
        '''Syntax: len(self)

        :return: The number of genes.'''
        return self.n_genes

    def __contains__(self, key):
        ''' Syntax: key in self

        Checks whether key is in self.

        :param key: May either be the gene name or the gene id'''
        return key in self._idx

    def remove_chromosome(self, chromosome):
        '''Deletes the chromosome from the transcriptome

        :param chromosome: Name of the chromosome to remove'''
        del self.data[chromosome]
        self.make_index()

    def _get_sample_idx(self, group_column='name'):
        'a dict with group names as keys and index lists as values'
        return self.infos['sample_table'].groupby(group_column).groups

    @property
    def sample_table(self):
        '''The sample table contains sample names, group information, long read coverage, as well as all other potentially
        relevant information on the samples.'''
        try:
            return self.infos['sample_table']
        except KeyError:
            return pd.DataFrame(columns=['name', 'file', 'group', 'nonchimeric_reads', 'chimeric_reads'], dtype='object')

    @property
    def samples(self) -> list:
        '''An ordered list of sample names.'''
        return list(self.sample_table.name)

    def groups(self, by='group') -> dict:
        '''Get sample groups as defined in columns of the sample table.

        :param by: A column name of the sample table that defines the grouping.
        :return: Dict with groupnames as keys and list of sample names as values.
        '''
        return dict(self.sample_table.groupby(by)['name'].apply(list))

    @property
    def n_transcripts(self) -> int:
        '''The total number of transcripts isoforms.'''
        if self.data is None:
            return 0
        return sum(g.n_transcripts for g in self)

    @property
    def n_genes(self) -> int:
        '''The total number of genes.'''
        if self.data is None:
            return 0
        return sum((len(t) for t in self.data.values()))

    @property
    def novel_genes(self) -> int:  # this is used for id assignment
        '''The total number of novel (not reference) genes.'''
        try:
            return self.infos['novel_counter']
        except KeyError:
            self.infos['novel_counter'] = 0
            return 0

    @property
    def chromosomes(self) -> list:
        '''The list of chromosome names.'''
        return list(self.data)

    def __str__(self):
        return '{} object with {} genes and {} transcripts'.format(type(self).__name__, self.n_genes, self.n_transcripts)

    def __repr__(self):
        return object.__repr__(self)

    def __iter__(self):
        return (gene for tree in self.data.values() for gene in tree)

    # IO: load new data from primary data files
    from ._transcriptome_io import add_sample_from_bam, remove_samples, add_short_read_coverage, remove_short_read_coverage, collapse_immune_genes

    # IO: output data as tables or other human readable format
    from ._transcriptome_io import gene_table, transcript_table, chimeric_table, write_gtf, export_alternative_splicing

    # filtering functionality and iterators
    from ._transcriptome_filter import add_qc_metrics, add_filter, remove_filter, iter_genes, iter_transcripts, iter_ref_transcripts

    # statistic: differential splicing, alternative_splicing_events
    from ._transcriptome_stats import altsplice_test, alternative_splicing_events, rarefaction

    # statistic: summary tables (can be used as input to plot_bar / plot_dist)
    from ._transcriptome_stats import altsplice_stats, filter_stats, transcript_length_hist, transcript_coverage_hist,\
        transcripts_per_gene_hist, exons_per_transcript_hist, downstream_a_hist, direct_repeat_hist
