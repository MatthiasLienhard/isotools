# Changelog

## [TODO]
* change: extend flanking exons for MISO/rMATS export 
* integrate reference and isoseq genes
    * Gene.get_splice_graph(force_recaculation=False, min_transcripts=None)
* implement some functionality of altsplice_test in find_splice_bubbles()
    * find_splice_bubbles(weights=g.coverage) ...
    * implement "novel" also add for PCA / find_splice_bubbles
    * implement filter for functions calling find_splice_bubbles by setting weights=0
* consider revising splice graph to common nomenclauture:
    * nodes/vertices are donor/acceptor sites
    * edges are exons/introns
    * pro:
        * standard
        * trivial to reconstruct transcripts
        * most functions should work (faster?)
        * can be extended easily
    * contra:
        * not trivial to find overlapping transcripts (not needed?)
        * breaks: is_same_exon,_check_exon,_check_junction, get_overlap



## [0.0.2] - 2020-03-05

* new feature: added decorators "experimental" and "deprecated" to mark unsave functions 
* change: in differential splicing changed the alternative fraction, to match the common PSI (% spliced in) definition
* change: narrowed definition of mutually exclusive exons: the alternatives now need to to feature exactly one ME exon and rejoin at node C
* change: for ME exons now the beginning of node C is returned as "end" of the splice bubble
* new feature: differential splicing result contains "novel", indicating that the the alternative is in the annotation 
* new feature: added alternative TSS/alternative PAS to the differential splicing test
* change: removed obsolete weights from splice graph and added strand
* change: unified parameters and column names of results of Transcriptome.find_splice_bubbles() and Transcriptome.altsplice_test()
* fix: add_short_read_coverage broken if short reads are already there. 


## [0.0.1] - 2020-02-25
* first shared version
* new feature: added option to export alternative splicing events for MISO and rMATS
* new feature: added changelog

