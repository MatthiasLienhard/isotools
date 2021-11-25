[![tests](https://github.com/MatthiasLienhard/isotools/actions/workflows/tests.yml/badge.svg)](https://github.com/MatthiasLienhard/isotools/actions?query=workflow%3Atests)
[![docs](https://readthedocs.org/projects/isotools/badge/?version=latest)](https://isotools.readthedocs.io/en/latest/)
[![PyPI](https://img.shields.io/pypi/v/isotools.svg)](https://pypi.org/project/isotools)
[![PyPIDownloadsTotal](https://pepy.tech/badge/isotools)](https://pepy.tech/project/isotools)
[![Licence: MIT](https://img.shields.io/badge/license-MIT-blue)](https://github.com/MatthiasLienhard/isotools/blob/master/LICENSE.txt)
# IsoTools 

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
This code block demonstrates the basic file import with IsoTools. 
It uses a small test data set contained in this repository, and should run within seconds. The paths are relative to the root of the repository.
For more comprehensive real world examples see the [tutorials](https://isotools.readthedocs.io/en/latest/tutorials.html "readthedocs").
```python
from  isotools import Transcriptome
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
# import the reference annotation
isoseq = Transcriptome.from_reference('tests/data/example.gff.gz')
# import the isoseq data
for sa in ('CTL', 'VPA'):
    isoseq.add_sample_from_bam(f'../tests/data/example_1_{sa}.bam', sample_name=sa, group=sa, platform='SequelII')
# save the imported file as pkl file (for faster import)
isoseq.add_qc_metrics('../tests/data/example.fa')
isoseq.save('../tests/data/example_1_isotools.pkl')
```

## Citation and feedback:
* If you run into any issues, please use the [github issues report feature](https://github.com/MatthiasLienhard/isotools/issues).
* For general feedback, please write me an email to [lienhard@molgen.mpg.de](mailto:lienhard@molgen.mpg.de).
* If you use isotools in your publication, please cite the following paper: IsoTools: IsoTools: a python toolbox for long-read transcriptome sequencing (in preparation)
