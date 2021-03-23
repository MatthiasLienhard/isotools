# Change Log

## [TODO] ideas and planed extensions or changes that are not yet implemented
* extend flanking exons for MISO/rMATS export 
* integrate reference and isoseq genes
    * Gene.get_splice_graph(force_recaculation=False, min_transcripts=None)
* keep track of actual (read) TSS and PAS within first/last exon
* implement some functionality of altsplice_test in find_splice_bubbles()
    * find_splice_bubbles(weights=g.coverage) ...
    * implement "novel" also add for PCA / find_splice_bubbles
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

