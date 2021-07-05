# isotools
IsoTools is a python module for Long Read Transcriptome Sequencing (LRTS) analysis.

Key features:
* Import of LRTS bam files (aligned full length transcripts).
* Import of reference annotation in gff3/gtf format.
* Computation of quality control metrics.
* Annotation and classification of novel transcripts with biologically motivated classification scheme.
* Definition of alternative splicing events based on segment graphs.
* Detection of differential alternative splicing between samples and groups of samples.
* Data visualization.

## documentation:
The documentation, including tutorials with real-world case studies and the complete API reference is available at [readthedocs](https://isotools.readthedocs.io/en/latest/ "documentation")

## installation:
Isotools is available from PyPI, and can be installed with the pip command:
```
python3 -m pip install isotools

```
Alternatively, to install from github, use the following command:

```
git clone https://github.com/MatthiasLienhard/isotools.git
cd isotools
python3 -m pip install .
```

## usage:
This code block demonstrates the basic file import with isoseq. For a more comprehensive real world example see the [tutorials](https://isotools.readthedocs.io/en/latest/tutorials.html "readthedocs").
```python
from  isotools.transcriptome import Transcriptome
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
```

## Citation and feedback:
* If you run into any issues, please use the [github issues report feature](https://github.com/MatthiasLienhard/isotools/issues).
* For general feedback, please write me an email to [lienhard@molgen.mpg.de](mailto:lienhard@molgen.mpg.de).
* If you use isotools in your publication, please cite the following paper: IsoTools: IsoTools: a python toolbox for long-read transcriptome sequencing (in preparation)
