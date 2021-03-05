from pysam import TabixFile, AlignmentFile, FastaFile
from tqdm import tqdm
import logging
logger=logging.getLogger('isotools')

default_gene_filter={'NOVEL_GENE':'not reference',
                    'EXPRESSED':'transcripts',
                    'CHIMERIC':'chimeric'}

default_ref_transcript_filter={
        'UNSPLICED':'len(exons)==1',
        'MULTIEXON':'len(exons)>1',
        'INTERNAL_PRIMING':'downstream_A_content and downstream_A_content>.5'}

default_transcript_filter={
        #'CLIPPED_ALIGNMENT':'clipping',
        'INTERNAL_PRIMING':'len(exons)==1 and downstream_A_content and downstream_A_content>.5', #more than 50% a
        'RTTS':'noncanonical_splicing and novel_splice_sites and any(2*i in novel_splice_sites[0] and 2*i+1 in novel_splice_sites[0] for i,_ in noncanonical_splicing)',
        'NONCANONICAL_SPLICING':'noncanonical_splicing',
        'NOVEL_TRANSCRIPT':'annotation is None or annotation[0]>0',
        'FRAGMENT':'fragments',
        'NOVEL':'not annotation or annotation[0]==4',
        'UNSPLICED':'len(exons)==1',
        'MULTIEXON':'len(exons)>1'}

annotation_vocabulary=['antisense', 'intergenic', 'genic genomic', 'novel exonic splice donor', 'novel exonic splice acceptor', 'novel PAS', 'readthrough fusion', 'novel exon', 'novel intronic splice donor', 'intron retention', 'novel intronic splice acceptor', 'exon skipping', 'novel combination', 'novel TSS', 'mono-exon', 'novel junction', "5' fragment", "3' fragment", 'intronic']
# filtering functions for the transcriptome class
def add_biases(self, genome_fn):
    'populates transcript["biases"] information, which can be used do create filters'
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
    'create filter flags which can be used by iter_transcripts'
    gene_attributes={k for g in self for k in g.data.keys() }
    tr_attributes={k for g in self for tr in g.transcripts for k in tr.keys() }
    ref_tr_attributes={k for g in self if g.is_annotated for tr in g.ref_transcripts for k in tr.keys() }
    tr_attributes.add('filter')
    ref_tr_attributes.add('filter')
    if  gene_filter is None:
        gene_filter=default_gene_filter
    if transcript_filter is None:
        transcript_filter=default_transcript_filter
    if ref_transcript_filter is None:
        ref_transcript_filter=default_ref_transcript_filter
    gene_ffun={label:_filter_function(gene_attributes, fun) for label,fun in gene_filter.items()}
    tr_ffun={label:_filter_function(tr_attributes, fun) for label,fun in transcript_filter.items()}
    reftr_ffun={label:_filter_function(ref_tr_attributes, fun) for label,fun in ref_transcript_filter.items()}
    for g in tqdm(self):
            g.add_filter(gene_ffun,tr_ffun,reftr_ffun)
    self.infos['filter']={'gene_filter':gene_filter, 'transcript_filter':transcript_filter, 'ref_transcript_filter':ref_transcript_filter}

def iter_genes(self, region=None,include=None, remove=None):
    'iterate over the genes of a region, optionally applying filters'
    if include or remove:
        assert 'filter' in self.infos, 'no filter flags found - run .add_filter() method first'
        assert not include or all(f in self.infos['filter']['gene_filter'] for f in include), 'not all filters to include found'
        assert not remove or all(f in self.infos['filter']['gene_filter'] for f in remove), 'not all filters to remove found'
    if region is None:
        genes=self
    elif region in self.data:
        genes=self.data[region] #provide chromosome
    else:
        try:
            chrom,start,end=region
        except ValueError:
            chrom, pos=region.split(':')
            start, end=[int(v) for v in pos.split('-')]
        except:
            raise ValueError('incorrect region {} - specify as string "chr:start-end" or tuple ("chr",start,end)'.format(region))
        finally:
            genes=self.data[chrom][start:end]
    for g in genes:
        if not include or any(f in include for f in g.data['filter']):
            if not remove or all(f not in remove for f in g.data['filter']):        
                yield g

def iter_transcripts(self,region=None,include=None, remove=None, min_coverage=None, max_coverage=None):
    'iterate over the transcripts of a region, optionally applying filters'   
    
    if include or remove:
        valid_filters=annotation_vocabulary
        if isinstance(include, str):
            include=[include]
        if isinstance(remove, str):
            remove=[remove]
        if 'filter' in self.infos:
            all_filter=self.infos['filter']
            valid_filters=valid_filters+list(all_filter['transcript_filter'])+list(all_filter['gene_filter'] )
        msg='did not find the following filter flags for {}: {}\nvalid filters are: {}'
        assert not include or all(f in valid_filters for f in include), msg.format( 
            'inclusion', ', '.join(f for f in include if f not in valid_filters), ', '.join(valid_filters) )
        assert not remove or all(f in valid_filters for f in remove), msg.format(
            'removal', ', '.join(f for f in remove if f not in valid_filters), ', '.join(valid_filters) )

    anno_include=[f for f in include if f in annotation_vocabulary] if include else []
    anno_remove=[f for f in remove if f in annotation_vocabulary] if remove else []
    g_include=[f for f in include if f in all_filter['gene_filter']] if include else []
    g_remove=[f for f in remove if f in all_filter['gene_filter']] if remove else []
    t_include=[f for f in include if f in all_filter['transcript_filter']] if include else []
    t_remove=[f for f in remove if f in all_filter['transcript_filter']] if remove else []
    
    for g in self.iter_genes(region, g_include, g_remove):
        for i,tr in g.filter_transcripts(t_include, t_remove, anno_include, anno_remove, min_coverage, max_coverage):
            yield g,i,tr

def iter_ref_transcripts(self,region=None,include=None, remove=None):
    'iterate over the transcripts of a region, optionally applying filters'   
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
    
