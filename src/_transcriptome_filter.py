from pysam import TabixFile, AlignmentFile, FastaFile
from tqdm import tqdm
import logging
logger=logging.getLogger('isotools')

DEFAULT_GENE_FILTER={'NOVEL_GENE':'not reference',
                    'EXPRESSED':'transcripts',
                    'CHIMERIC':'chimeric'}
'''Default definitions for gene filter, as used in iosotools.Transcriptome.add_filters().'''

DEFAULT_REF_TRANSCRIPT_FILTER={
        'UNSPLICED':'len(exons)==1',
        'MULTIEXON':'len(exons)>1',
        'INTERNAL_PRIMING':'downstream_A_content and downstream_A_content>.5'}
'''Default definitions for reference transcript filter, as used in iosotools.Transcriptome.add_filters().'''

#:Default definitions for transcript filter, as used in iosotools.Transcriptome.add_filters()
DEFAULT_TRANSCRIPT_FILTER={
        #'CLIPPED_ALIGNMENT':'clipping',
        'INTERNAL_PRIMING':'len(exons)==1 and downstream_A_content and downstream_A_content>.5', #more than 50% a
        'RTTS':'noncanonical_splicing and novel_splice_sites and any(2*i in novel_splice_sites[0] and 2*i+1 in novel_splice_sites[0] for i,_ in noncanonical_splicing)',
        'NONCANONICAL_SPLICING':'noncanonical_splicing',
        'NOVEL_TRANSCRIPT':'annotation is None or annotation[0]>0',
        'FRAGMENT':'fragments and any("novel exonic " in a or "fragment" in a for a in annotation[1])' ,
        'NOVEL':'not annotation or annotation[0]==4',
        'UNSPLICED':'len(exons)==1',
        'MULTIEXON':'len(exons)>1'}


ANNOTATION_VOCABULARY=['antisense', 'intergenic', 'genic genomic', 'novel exonic PAS','novel intronic PAS', 'readthrough fusion', 
'novel exon', "novel 3' splice site", 'intron retention', "novel 5' splice site", 'exon skipping', 'novel combination', 
'novel intronic TSS','novel exonic TSS', 'mono-exon', 'novel junction', "5' fragment", "3' fragment", 'intronic']
'''Controlled vocabulary for filtering by novel alternative splicing.'''

# filtering functions for the transcriptome class
def add_qc_metrics(self, genome_fn):
    ''' Retrieves QC metrics for the transcripts. 

    Calling this function populates transcript["biases"] information, which can be used do create filters. 
    In particular, the direct repeat length, the downstream adenosine content and information about noncanonical splice sites are fetched.
    Additionaly genes are scanned for transcripts that are fully contained in other transcripts. 
    
    :param geneome_fn: path to the genome in fastA format.'''
    with FastaFile(genome_fn) as genome_fh:
        missing_chr=set(self.chromosomes)-set(genome_fh.references)
        if missing_chr:
            missing_genes=sum(len(self.data[mc]) for mc in missing_chr)
            logger.warning(f'{len(missing_chr)} contigs are not contained in genome, affecting {missing_genes} genes. Some metrics cannot be computed: {missing_chr}')

        for g in tqdm(self):                
            g.add_fragments()
            if g.chrom in genome_fh.references:
                g.add_direct_repeat_len(genome_fh) 
                g.add_noncanonical_splicing(genome_fh)
                g.add_threeprime_a_content(genome_fh)
                
    self.infos['biases']=True # flag to check that the function was called

def add_filter(self, gene_filter=None,transcript_filter=None, ref_transcript_filter=None):
    '''Defines and assigns filter flags, which can be used by iter_transcripts.
    
    Filters are defined as dict, where the key is a filter identifier, and the value is an expression, 
    which gets evaluated on the gene/transcript. For examples, see the default filter definitions 
    isotools.DEFAULT_GENE_FILTER, isotools.DEFAULT_TRANSCRIPT_FILTER and isotools.DEFAULT_REF_TRANSCRIPT_FILTER.

    :param gene_filter: dict of gene filters. If omitted the default gene filters apply.
    :param transcript_filter: dict of gene filters. If omitted the default reference filters apply.
    :param ref_transcript_filter: dict of gene filters. If omitted the default transcript filters apply.
    '''
    gene_attributes={k for g in self for k in g.data.keys() if k.isidentifier()}
    tr_attributes={k for g in self for tr in g.transcripts for k in tr.keys() if k.isidentifier()}
    ref_tr_attributes={k for g in self if g.is_annotated for tr in g.ref_transcripts for k in tr.keys() if k.isidentifier()}
    tr_attributes.add('filter')
    ref_tr_attributes.add('filter')
    if  gene_filter is None:
        gene_filter=DEFAULT_GENE_FILTER
    if transcript_filter is None:
        transcript_filter=DEFAULT_TRANSCRIPT_FILTER
    if ref_transcript_filter is None:
        ref_transcript_filter=DEFAULT_REF_TRANSCRIPT_FILTER
    gene_ffun={label:_filter_function(gene_attributes, fun) for label,fun in gene_filter.items()}
    tr_ffun={label:_filter_function(tr_attributes, fun) for label,fun in transcript_filter.items()}
    reftr_ffun={label:_filter_function(ref_tr_attributes, fun) for label,fun in ref_transcript_filter.items()}
    for g in tqdm(self):
            g.add_filter(gene_ffun,tr_ffun,reftr_ffun)
    self.infos['filter']={'gene_filter':gene_filter, 'transcript_filter':transcript_filter, 'ref_transcript_filter':ref_transcript_filter}

def iter_genes(self, region=None,include=None, remove=None):
    '''Iterates over the genes of a region, optionally applying filters.
    
    :param region: The region to be considered. Either a string "chr:start-end", or a tuple (chr,start,end). Start and end is optional. 
    :param include: If provided, only genes featuring at least one of these flags are considered.
    :param remove: If provided, genes featuring one of at least these flags are ignored.'''

    if include or remove:
        assert 'filter' in self.infos, 'no filter flags found - run .add_filter() method first'
        assert not include or all(f in self.infos['filter']['gene_filter'] for f in include), 'not all filters to include found'
        assert not remove or all(f in self.infos['filter']['gene_filter'] for f in remove), 'not all filters to remove found'
    if region is None:
        genes=self
    elif isinstance(region, str):
        if region in self.data:
            genes=self.data[region] #provide chromosome
        else:
            try:
                chrom, pos=region.split(':')
                if chrom in self.data:
                    start, end=[int(v) for v in pos.split('-')]
                else:
                    raise ValueError('specified chromosome {} not found'.format(chrom))
            except:
                raise ValueError('incorrect region {} - specify as string "chr:start-end" or tuple ("chr",start,end)'.format(region))
            if chrom in self.data:
                genes=self.data[chrom][int(start):int(end)]
            else:
                raise ValueError('specified chromosome {} not found'.format(chrom))
    elif isinstance(region,tuple):
        chrom,start,end=region
        if chrom in self.data:
            genes=self.data[chrom][int(start):int(end)]
        else:
            raise ValueError('specified chromosome {} not found'.format(chrom))
    for g in genes:
        if not include or any(f in include for f in g.data['filter']):
            if not remove or all(f not in remove for f in g.data['filter']):        
                yield g

def iter_transcripts(self,region=None,include=None, remove=None, min_coverage=None, max_coverage=None):
    '''Iterates over the transcripts of a region, optionally applying filters.
    :param region: The region to be considered. Either a string "chr:start-end", or a tuple (chr,start,end). Start and end is optional. 
    :param include: If provided, only transcripts featuring at least one of these flags are considered.
    :param remove: If provided, transcripts featuring one of at least these flags are ignored.
    :param min_coverage: The minimum coverage threshold. Transcripts with less reads are ignored. 
    :param max_coverage: The maximum coverage threshold. Transcripts with more reads are ignored. '''
    if include or remove:
        valid_filters=ANNOTATION_VOCABULARY
        if isinstance(include, str):
            include=[include]
        if isinstance(remove, str):
            remove=[remove]
        if 'filter' in self.infos:
            all_filter=self.infos['filter']
        else:
            all_filter={'transcript_filter':{}, 'gene_filter':{}}
        valid_filters=valid_filters+list(all_filter['transcript_filter'])+list(all_filter['gene_filter'] )
        msg='did not find the following filter flags for {}: {}\nvalid filters are: {}'
        assert not include or all(f in valid_filters for f in include), msg.format( 
            'inclusion', ', '.join(f for f in include if f not in valid_filters), ', '.join(valid_filters) )
        assert not remove or all(f in valid_filters for f in remove), msg.format(
            'removal', ', '.join(f for f in remove if f not in valid_filters), ', '.join(valid_filters) )

    anno_include=[f for f in include if f in ANNOTATION_VOCABULARY] if include else []
    anno_remove=[f for f in remove if f in ANNOTATION_VOCABULARY] if remove else []
    g_include=[f for f in include if f in all_filter['gene_filter']] if include else []
    g_remove=[f for f in remove if f in all_filter['gene_filter']] if remove else []
    t_include=[f for f in include if f in all_filter['transcript_filter']] if include else []
    t_remove=[f for f in remove if f in all_filter['transcript_filter']] if remove else []
    
    for g in self.iter_genes(region, g_include, g_remove):
        for i,tr in g.filter_transcripts(t_include, t_remove, anno_include, anno_remove, min_coverage, max_coverage):
            yield g,i,tr

def iter_ref_transcripts(self,region=None,include=None, remove=None):
    '''Iterates over the referemce transcripts of a region, optionally applying filters.
    
    :param region: The region to be considered. Either a string "chr:start-end", or a tuple (chr,start,end). Start and end is optional. 
    :param include: If provided, only genes featuring at least one of these flags are considered.
    :param remove: If provided, genes featuring one of at least these flags are ignored.'''
    if include or remove:
        assert 'filter' in self.infos, 'no filter flags found - run .add_filter() method first'
        all_filter=self.infos['filter']
    g_include=[f for f in include if f in all_filter['gene_filter']] if include else []
    g_remove=[f for f in remove if f in all_filter['gene_filter']] if remove else []
    t_include=[f for f in include if f not in all_filter['gene_filter']] if include else []
    t_remove=[f for f in remove if f not in all_filter['gene_filter']] if remove else []
    assert all(f in all_filter['ref_transcript_filter'] for f in t_include), 'not all reference filters to include found'
    assert all(f in all_filter['ref_transcript_filter'] for f in t_remove), 'not all  reference filters to remove found'
    for g in self.iter_genes(region, g_include, g_remove):
        if g.is_annotated:
            for i,tr in g.filter_ref_transcripts(t_include, t_remove):
                yield g,i,tr


def _filter_function(argnames, expression):
    'converts a string e.g. "all x[0]/x[1]>3" into a function'
    return eval (f'lambda {",".join(arg+"=None" for arg in argnames)}: bool({expression})\n',{},{})
    
