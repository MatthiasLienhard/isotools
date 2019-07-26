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


class Transcriptome(dict):
    def __init__(self, bam_fn=None, gtf_fn=None,**kwargs):
        if bam_fn is not None and gtf_fn is not None:
            raise ValueError('specify either bam_fn or gtf_fn')