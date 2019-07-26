# isotools
python module for isoseq postprocessing
* import of bam files (aligned isoseq transcripts)
* collapsing isoseq transcripts
    * truncation
    * fuzzy junctions (todo)
* import of gff files (e.g. for reference annotation, but also isoseq in gtf format)
* comparison to reference annotation (e.g. refseq)
    * for gene names
    * to detect alternative splicing wrt referernce

## TODO:
* class implementation:
    * transcriptome.Transcriptome()
* collapse isoforms:
    * fuzzy junctions
    * thresholds for truncation
    * consider misalignments???

