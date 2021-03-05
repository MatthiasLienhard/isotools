# Changelog

## [TODO]
* change: extend flanking exons for MISO/rMATS export 


## [Unreleased]

* new feature: added decorators "experimental" and "deprecated" to mark unsave functions 
* change: in differential splicing changed the alternative fraction, to match the common PSI (% spliced in) definition
* change: narrowed definition of mutually exclusive exons: the alternatives now need to to feature exactly one ME exon and rejoin at node C
* change: for ME exons now the beginning of node C is returned as "end" of the splice bubble
* new feature: differential splicing result contains "novel", indicating that the the alternative is in the annotation (todo: also add for PCA / find_splice_bubbles)
* new feature: added alternative TSS/alternative PAS to the differential splicing test
* change: removed obsolete weights from splice graph and added strand
* change: unified parameters and column names of results of Transcriptome.find_splice_bubbles() and Transcriptome.altsplice_test()
* fix: add_short_read_coverage broken if short reads are already there. 


## [0.0.1] - 2020-02-25
* first shared version
* new feature: added option to export alternative splicing events for MISO and rMATS
* new feature: added changelog

