# Change Log

## [TODO] ideas and planed extensions or changes that are not yet implemented
* run_isotools console script is totally outdated and broken...
* optimize add_biases for run after new samples have been added - should not recompute everything
* extend flanking exons for MISO/rMATS export (not really needed, works fine as is)
* consider integrating reference and isoseq genes one in segment graph
    * Gene.get_segment_graph(force_recaculation=False, min_transcripts=None)
    * would be major change and potentially break things
* keep track of actual (read) TSS and PAS within first/last exon
* implement some functionality of altsplice_test in find_splice_bubbles()
    * find_splice_bubbles(weights=g.coverage) ...
    * implement filter for functions calling find_splice_bubbles by setting weights=0
* consider using common splice graph e.g. to import bam:
    * nodes/vertices are donor/acceptor sites
    * edges are exons/introns
    * pro:
        * commonly used
        * trivial to reconstruct transcripts
        * most functions should work (faster?)
        * can be extended easily
    * con:
        * segment graph still needed for bubble definition - two graphs stored 

## [0.1.2]

* New: added function remove_short_read_coverage
* New: added some missing documentation for gene plots
* Fix: fixed bug in novel transcript class definition, affecting last exons
* New: Distinguish novel exonic TSS (NIC) and novel intronic TSS (NNC)
* New: Do not distinguish intronic/exonic novel splice sites. Report distance to shortest splice site of same type.
* Fix: Sashimi plots ignored mono exons


## [0.1.1] - 2020-04-12

* Fix: fixed bug in TSS/PAS events affecting start/end positions and known flag.
* Change: refactored Transcriptome.find_splice_bubbles() to Transcriptome.alternative_splicing_events()
* Change: refactored SegmentGraph.find_alternative_starts() to SegmentGraph.find_start_end_events()

## [0.1.0] - 2020-03-24

* added documentation
* moved examples in documentation

## [0.0.2] - 2020-03-22
* Change: refactored SpliceGraph to SegmentGraph to better comply with common terms in literature
* New: added a basic implementation of an actual SpliceGraph (as commonly defined in literature) 
    * based on sorted dict
    * not used so far, but maybe useful in importing the long read bam files since it can be extended easily
* New: added decorators "experimental" and "deprecated" to mark unsafe functions 
* Change: in differential splicing changed the alternative fraction, to match the common PSI (% spliced in) definition
* Change: narrowed definition of mutually exclusive exons: the alternatives now need to to feature exactly one ME exon and rejoin at node C
* Change: for ME exons now the beginning of node C is returned as "end" of the splice bubble
* New: differential splicing result contains "novel", indicating that the the alternative is in the annotation 
* New: added alternative TSS/alternative PAS to the differential splicing test
* Change: removed obsolete weights from splice graph and added strand
* Change: unified parameters and column names of results of Transcriptome.find_splice_bubbles() and Transcriptome.altsplice_test()
* Fix: add_short_read_coverage broken if short reads are already there. 


## [0.0.1] - 2020-02-25
* first shared version
* New: added option to export alternative splicing events for MISO and rMATS
* New: added change log

