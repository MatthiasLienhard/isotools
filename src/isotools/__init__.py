
import sys
if sys.version_info >= (3, 8):
    from importlib.metadata import distribution
else:
    from importlib_metadata import distribution  # py3.7

__version__ = distribution('isotools').version

from ._transcriptome_filter import DEFAULT_GENE_FILTER, DEFAULT_TRANSCRIPT_FILTER, DEFAULT_REF_TRANSCRIPT_FILTER, ANNOTATION_VOCABULARY
from .splice_graph import SegmentGraph, SegGraphNode
from .transcriptome import Transcriptome
from .gene import Gene


__all__ = ['Transcriptome', 'Gene', 'SegmentGraph', 'SegGraphNode',
           'DEFAULT_GENE_FILTER', 'DEFAULT_TRANSCRIPT_FILTER', 'DEFAULT_REF_TRANSCRIPT_FILTER', 'ANNOTATION_VOCABULARY']
