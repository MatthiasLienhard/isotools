Getting Started
===============
Python module for Long Read Transcriptome Sequencing (LRTS) quality control and analysis.

Key Features:

* Import of LRTS bam files (aligned full length transcripts).
* Import of reference annotation in gff3/gtf format.
* Computation quality control metrics.
* Annotation and classification of novel transcripts with biologically motivated classification scheme.
* Definition of alternative splicing events based on segment graphs.
* Detection of differential alternative splicing between samples and groups of samples. 
* Data visualization. 

Installation
------------
To install the package, download the latest version from github and install with pip:

.. code-block:: bash

    git clone https://github.molgen.mpg.de/lienhard/isotools.git
    cd isotools
    pip install .

Usage
-----
This code block demonstrates the basic file import with isoseq. For a more comprehensive real world example see the tutorial. 

.. code-block:: python

    from  isotools import Transcriptome
    import logging
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    isoseq=Transcriptome.from_reference('reference_file.gff3.gz')
    isoseq_bam_fn={'sample1':'isoseq_fn_s1.bam', 'sample2':'isoseq_fn_s2.bam'}
    groups={'sample1':'control', 'sample2':'treatment'}
    for sa,bam in isoseq_bam_fn.items():
        isoseq.add_sample_from_bam(bam, sample_name=sa, group=groups[sa]) 
    isoseq.add_biases('genome.fa')
    isoseq.make_index()
    isoseq.add_filter()
    isoseq.save('example_isotools.pkl')

Citation
--------
If you use isotools in your publication, please cite the following paper:

IsoTools: the python toolbox for characterization,expression quantification, 
and differential analysis of isoforms from long read sequencing experiments (in preparation)