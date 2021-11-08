
try:
    from importlib.metadata import distribution
except ModuleNotFoundError:
    from importlib_metadata import distribution  # py3.7
__version__ = distribution('isotools').version
from .gene import Gene
from .transcriptome import Transcriptome
from .splice_graph import SegmentGraph, SegGraphNode

from ._transcriptome_filter import DEFAULT_GENE_FILTER, DEFAULT_TRANSCRIPT_FILTER, DEFAULT_REF_TRANSCRIPT_FILTER, ANNOTATION_VOCABULARY


__all__ = ['Transcriptome', 'Gene', 'SegmentGraph', 'SegGraphNode',
           'DEFAULT_GENE_FILTER', 'DEFAULT_TRANSCRIPT_FILTER', 'DEFAULT_REF_TRANSCRIPT_FILTER', 'ANNOTATION_VOCABULARY']
