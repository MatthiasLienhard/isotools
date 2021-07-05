# Change Log

## TODO: ideas, issues and planed extensions or changes that are not yet implemented
* avoid the need for add_filters - construct and evaluate lambdas during filtering
* run_isotools console script is totally outdated and broken
* tests are outdated/broken/not automated
* optimize add_qc_metrics for run after new samples have been added - should not recompute everything


## [0.2]
* restructure to meet PyPI recommendations

## [0.1.5]
* New feature: restrict tests on provided splice_types
* New feature: provide position to find given alternative splicing events



## [0.1.4]
* Fix: Issue with noncanonical splicing detection introduced in 0.1.3
* Fix: crash with secondary alignments in bam files during import.
* New feature: Report and skip if alignment outside chromosome (uLTRA issue)
* Fix: import of chimeric reads (secondary alignments have no SA tag)
* Fix: Transcripts per sample in sample table: During import count only used transcripts, do not count chimeric transcripts twice. 
* Change: sample_table reports chimeric_reads and nonchimeric_reads (instead of total_reads)
* Change: import of long read bam is more verbose in info mode
* Fix: Bug: import of chained chimeric alignments overwrites read coverage when merging to existing transcript
* Fix: remove_samples actually removes the samples from the sample_table
* Change: refactored add_biases to add_qc_metrics
* fix: property of transcripts included {sample_name:0}
* save the TSS and PAS positions
* New: use_satag parameter for add_sample_from_bam 
* Change: use median TSS/PAS (of all reads with same splice pattern) as transcript start/end (e.g. exons[0][0]/exons[-1][1])
* Fix: Novel exon skipping annotation now finds all exonic regions that are skipped.
* change: Default filter of FRAGMENTS now only tags reads that do not use a reference TSS or PAS
## [0.1.3]
* Fix: improved performance of noncanonical splicing detection by avoiding redundant lookups. 


## [0.1.2] - 2020-05-03

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

