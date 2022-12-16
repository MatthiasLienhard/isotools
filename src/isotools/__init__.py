"""
IsoTools: Python package for long read transcriptome sequencing analysis.

.. data:: DEFAULT_GENE_FILTER

    Default definitions for gene filter, as used in iosotools.Transcriptome.add_filters().

.. data:: DEFAULT_TRANSCRIPT_FILTER

    Default definitions for transcript filter, as used in iosotools.Transcriptome.add_filters().

.. data:: DEFAULT_REF_TRANSCRIPT_FILTER

    Default definitions for reference transcript filter, as used in iosotools.Transcriptome.add_filters().

.. data:: ANNOTATION_VOCABULARY

    Controlled vocabulary for filtering by annotation.

.. data:: SPLICE_CATEGORY

    Controlled vocabulary for filtering by novel alternative splicing.
"""


try:
    from importlib.metadata import distribution
except ModuleNotFoundError:
    from importlib_metadata import distribution  # py3.7
__version__ = distribution('isotools').version
from .gene import Gene
from .transcriptome import Transcriptome
from .splice_graph import SegmentGraph, SegGraphNode
from ._transcriptome_stats import estimate_tpm_threshold

from ._transcriptome_filter import DEFAULT_GENE_FILTER, DEFAULT_TRANSCRIPT_FILTER, DEFAULT_REF_TRANSCRIPT_FILTER, ANNOTATION_VOCABULARY


__all__ = ['Transcriptome', 'Gene', 'SegmentGraph', 'SegGraphNode', 'estimate_tpm_threshold',
           'DEFAULT_GENE_FILTER', 'DEFAULT_TRANSCRIPT_FILTER', 'DEFAULT_REF_TRANSCRIPT_FILTER', 'ANNOTATION_VOCABULARY']
