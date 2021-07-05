# isotools
python module for isoseq postprocessing


## documentation:
http://medips.molgen.mpg.de/isoseq/


## installation:
Isotools is available from PyPI, and can be installed with the pip command:
```
python3 -m pip install isotools

```
To install from github:

```
git clone https://github.molgen.mpg.de/lienhard/isotools.git
cd isotools
python3 -m pip install .
```

## usage:
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
