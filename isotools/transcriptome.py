#!/usr/bin/python3
import os
import itertools
import warnings
from itertools import combinations

import argparse
from pysam import TabixFile, AlignmentFile
import pandas as pd
import numpy as np
#from Bio import SeqIO, Seq, SeqRecord
import matplotlib.pyplot as plt
from intervaltree import IntervalTree, Interval
from tqdm import tqdm
import isotools.io
import isotools.splice_graph
import isotools.utils

class Transcriptome():
    def __init__(self, fn=None, **kwargs):
        self.data=None
        if fn is not None:
            self.import_file(fn, **kwargs)
    
    def import_file(self, fn, **kwargs):
        self.data=isotools.io.import_transcripts(fn, **kwargs)

    def __len__(self):
        return self.n_transcripts
    
    @property
    def n_transcripts(self):
        if self.data==None:
            return 0
        return sum((len(g.data['transcripts']) for t in self.data.values() for g in t))

    @property
    def n_genes(self):
        if self.data==None:
            return 0
        return sum((len(t) for t in self.data.values()))
    
    def __str__(self):
        return '<{} object with {} genes and {} transcripts'.format(type(self).__name__, self.n_genes, self.n_transcripts)


    
def collapse_transcripts(self,  fuzzy_junction=5, repair_5truncation=10000, repair_3trunkation=100, rename=True):
    with tqdm(total=self.n_genes, unit='genes') as pbar:     
        for chrom, tree in self.data.items():
            for gene in tree:
                isotools.collapse.collapse_transcript_of_gene(gene, fuzzy_junction, repair_5truncation, repair_3trunkation, rename)
                pbar.update(1)

       
    def add_splice_graphs(self, force=False):
        with tqdm(total=self.n_genes, unit='genes') as pbar:     
            for chrom, tree in self.data.items():
                pbar.set_postfix(chr=chrom)
                for gene in tree:
                    pbar.update(1)
                    if not force and 'splice_graph' in gene.data:
                        continue
                    gene.data['splice_graph']=isotools.splice_graph.get_splice_graph(gene)


    
    def add_reference_support(self, ref_genes):
        n = sum([len(gene.data['transcripts'])
                    for tree in self.data.values() for gene in tree])
        with tqdm(total=n, unit='transcripts') as pbar:
            for chrom, tree in self.data.items():
                pbar.set_postfix(chr=chrom)
                for gene in tree:
                    for tr in gene.data['transcripts'].values():
                        tr['support']=isotools.utils.get_support(tr['exons'], ref_genes, chrom=chrom,is_reverse=gene.data['strand'] == '-')
                        pbar.update(1)

