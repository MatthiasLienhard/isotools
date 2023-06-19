# Change Log

## TODO: ideas, issues and planed extensions or changes that are not yet implemented
* optimize add_qc_metrics for run after new samples have been added - should not recompute everything
* planned new feature: during import of long reads, (optionally) correct for short exon alignment issues. 
* separate new read import and classification of isoforms.

## [0.3.4]
* fixing #8: AssertationError when unifying TSS/PAS between transcript
* improved domain plots: ORF start and end do not appear like exon exon boundaries. 
* API change: separated ORF prediction from QC metrics calculation. 
* new feature: count number of upstream start codons in Gene.add_orfs() (called by default when adding QC metrics to transcriptome)
* new feature: calculate Fickett testcode and hexamer score for longest ORFs, to separate coding and noncoding genes. 


## [0.3.3]
* fixed bug in filter_ref_transcripts with no query
* export gtf with long read transcripts as well uncovered as reference transcripts
* fix warning in plot_diff_results
* changed export to rMATS: events report complete flanking exons of top covered isoform
* fixed bug with transcript_id_col parameter in add_sample_from_csv
* fixed handling of interpro protein domains
* improved documentation: syntax highlighting, code style, additional explanations on filtering

## [0.3.2]
* restructured tutorials
* new feature: add domains to differential splicing result tables.
* new feature: min_coverage and max_coverage for iter_genes function.

## [0.3.1]
* new feature: add protein domains from 3 different sources and depict them with Gene.plot_domains()
* new feature: restrict gene and transcript iterators on list of genes of interest
* new feature: filter_transcripts function for genes
* changed SUBSTANTIAL filter to 1% of the genes total (was 5%)
* coordination test:
    * changed argument and column naming, to make it consistent with other test results
    * added conditional delta PSI effect size measure
    * order of events is now according to gene strand: A upstream of B

## [0.3.0]
* new feature: find longest ORF and infer NMD of lr transcripts (and annotation)
* new feature: allow for several TSS/PAS per intron chain and unify them across intron chains
* changed default parameter of filter_query in run_isotools script to "FSM or not (INTERNAL_PRIMING or RTTS)"

## [0.2.11.1]
* bugfix: KeyError during transcriptome reconstruction in _add_chimeric. 
* bugfix: default colors in plot_diff_results.

## [0.2.11]
* added function to import samples from csv/gtf to import transcriptome reconstruction / quantification from other tools.
* dropped requirement for gtf files to be tabix indexed.


## [0.2.10]
* fixed get_overlap - important for correct assignment of mono exonic genes to reference
* added parameter to control for minimal mapping quality in add_sample_from_bam. This allows for filtering out ambiguous reads, which have mapping quality of 0
* fixed plot_diff_result (Key error due to incorrect parsing of group names)
* New function estimate_tpm_threshold, to estimate the minimal abundance level of observable transcripts, given a sequencing depth. 
* New function coordination_test, to test coordination of splicing events within a gene. 
* Optional log or linear scale for the coverage axis in sashimi plots.

## [0.2.9]
* added DIE test
* adjusted classification of novel exonic TSS/PAS to ISM
* improved assignment of reference genes in case of equal number of matching splice sites to several reference genes. 
* added parameter to control for minimal exonic overlap to reference genes in add_sample_from_bam.
* changed computation of direct repeats. Added wobble and max_mm parameters.
* exposed parameters to end user in the add_qc_metrics function. 
* added options for additional fields in gtf output
* improved options for graphical output with the command line script
* fixed plot_bar default color scheme

## [0.2.8]
* fix: version information lost when pickeling reference.
* fix missing gene name
* added pt_size parameter to plot_embedding and plot_diff_results function
* added colors parameter to plotting functions
* various fixes of command line script run_isotools.py


## [0.2.7]
* added command line script run_isotools.py
* added test data for unit tests 


## [0.2.6]
* Added unit tests
* Fixed bug in novel splicing subcategory assignment
* new feature: rarefaction analysis
* Changed filtering: expressions get evaluated during iteration
    * Predefined filters are added automatically
    * Add / remove filters one by one
    * added optional progress bar to iter_genes/transcripts

## [0.2.5]
* New feature: distinguish noncanonical and canonical novel splice sites for direct repeat hist
* New feature: option to drop partially aligned reads with the min_align_fraction parameter in add_sample_from_bam

## [0.2.4]
* New feature: added option to save read names during bam import
* new feature: gzip compressed gtf output

## [0.2.3]
* Changed assignment of transcripts to genes if no splice sites match.
* Fix: more flexible import of reference files, gene name not required (but id is), introducing "infer_genes" from exon entries of gtf files.
* New function: Transcriptome.remove_filter(filter=[tags])

## [0.2.2]
* Fix: export to gtf with filter features

## [0.2.1]
* Fix: import reference from gtf file
* New feature: Import multiple samples from single bam tagged by barcode (e.g. from single cell data)
* Fix: issue with zero base exons after shifting fuzzy junctions


## [0.2.0]
* restructure to meet PyPI recommendations
* New feature: isoseq.altsplice_test accepts more than 2 groups, and computes ML parameters for all groups

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

