#!/usr/bin/env python
# coding: utf-8

# # Analysis workflow of the pacbio 2019 alzheimer Sequel2 dataset
# 
# ## Preparation
# 1) prepare the working directory and download the reference and data
# ``` sh
#     cd /my/working/directory
#     # create some subdirectories
#     mkdir -p reference isoseq/flnc isoseq/aligned isoseq/pickle tables plots
#     
#     # download a reference genome (806 MB)
#     genome_link='ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz'
#     wget -P reference -O GRCh38_genome_primary_assembly_gencode_v29.fa.gz ${genome_link} 
#     gunzip reference/GRCh38_genome_primary_assembly_gencode_v29.fa.gz
#     
#     # download gencode reference annotation (46.2 MB)
#     annotation_link='ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gff3.gz'
#     wget -P reference -O gencode_v29_annotation.gff3.gz ${annotation_link} 
#     
#     # download the isoseq flnc data (4.1 GB)
#     isoseq_link='https://downloads.pacbcloud.com/public/dataset/Alzheimer2019_IsoSeq/FullLengthReads/flnc.bam'
#     wget -P isoseq/flnc -O alzheimer2019_isoseq_flnc.bam ${isoseq_link} 
# ```
# 
# 2) align the isoseq data to the reference genome using minimap2.
# If the pacbio isoseq3 workflow is [installed](https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki/Tutorial:-Installing-and-Running-Iso-Seq-3-using-Conda) you can use the pacbio version of minimap2 as follows:
# 
# ``` sh
#     #activate the isoseq3 environement (assuming it is called pacbio)
#     conda activate pacbio
#     n_threads=8
#     sample='alzheimer2019_isoseq'
#     ref='reference/GRCh38_genome_primary_assembly_gencode_v29.fa'
#     pbmm2 align ${ref} isoseq/flnc/${sample}_flnc.bam isoseq/aligned/${sample}_aligned.sorted.bam --preset ISOSEQ --sort -j $n_threads 
# ```
# 
# 

# ## Isotools workflow

# In[1]:


from  isotools.transcriptome import Transcriptome
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# In[2]:


sample='alzheimer2019_isoseq'
bam=f'isoseq/aligned/{sample}_aligned.sorted.bam'
genome='reference/GRCh38_genome_primary_assembly_gencode_v29.fa'
anno='reference/gencode_v29_annotation'


# In[ ]:


try:
    #try to restore previously prepared data
    isoseq=Transcriptome(f'isoseq/pickle/{sample}_isotools.pkl')
except FileNotFoundError:
    #import the reference 
    isoseq=Transcriptome.from_reference(anno+'.gff3.gz') 
    # save the reference, so it it can be restored for other analysis
    isoseq.save_reference(anno+'.isotools.pkl') 
    # import the long read data and compare to reference (only one sample for this dataset)
    isoseq.add_sample_from_bam(bam, sample_name='alzheimer_1', group='alzheimer') 
    # compute QC metrics for all transcripts
    isoseq.add_biases(genome) #3:41 min
    # update the index for fast gene access by name/id
    isoseq.make_index()
    # save the current data
    isoseq.save(f'isoseq/pickle/{sample}_isotools.pkl')


# In[ ]:





# In[ ]:




