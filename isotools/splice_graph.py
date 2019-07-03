#!/usr/bin/python3

from pysam import TabixFile
import pysam
import pandas as pd
import numpy as np
from functools import reduce 
import operator
from Bio import SeqIO, Seq, SeqRecord
import csv
import argparse
import os
import itertools
import matplotlib.pyplot as plt
from intervaltree import IntervalTree, Interval
import warnings
from tqdm import tqdm
import argparse

class SpliceGraph:
    def __init__(self, refseq_fn="/project/42/references/refseq/RefSeq_GRCh38_20181116_sorted.gff.gz"):
        self.genes=dict() #one IntervalTree per chromosome, nodes contain gene position and infos
        print('read gff file {}'.format(refseq_fn))        
        gff = TabixFile(refseq_fn)        
        self.import_gtf(gff) #populate self.genes
        self.add_splice_graphs() #add the splice graph to self.genes
    
    @property
    def chrom(self):
        return(list(self.genes.keys()))

    @staticmethod
    def get_chrom_dict(gff):
        #fetch chromosome ids
        chrom={}
        for c in gff.contigs:
            #print ("---"+c)
            for line in gff.fetch(c,1,2):
                if line[1] =="C":
                    ls=line.split(sep="\t")
                    if ls[2] == "region":
                        info=dict([pair.split("=") for pair in ls[8].split(";")])
                        if "chromosome" in info.keys():
                            chrom[ls[0]]=info["chromosome"]
                            #print(line)
        return(chrom)     

    def import_gtf(self,gff):
        # returns a dict interval trees for the genes, each containing the splice graph
        chrom_ids=self.get_chrom_dict(gff)
        exons=dict() # transcript id -> exons
        transcripts=dict() # gene_id -> transcripts        
        skipped=set()
        for line in gff.fetch( ):#parser=pysam.asGTF()): parser is strange and not helpfull 
            ls=line.split(sep="\t")
            if ls[0] not in chrom_ids:
                #warnings.warn('unknown chromosome id :'+ls[0])
                #this happens at all NW/NT scaffolds
                continue
            chrom=chrom_ids[ls[0]]
            if chrom not in self.genes:
                self.genes[chrom]=IntervalTree()
            info=dict([pair.split("=",1) for pair in ls[8].split(";")])
            start, end=[int(i) for i in ls[3:5]]
            end+=1
            if 'Name' not in info:
                try:
                    info['Name']=info[{'standard_name','gene','product'}.intersection(info).pop()]
                except KeyError:
                    pass                
            if ls[2] == "exon" :
                #print(line)
                if "Parent" in info.keys():
                    gtf_id=info['Parent']
                    try:
                        exons[gtf_id].append((start, end))
                    except KeyError:
                        exons[gtf_id]=[(start, end)]
                else: #should not happen if GTF is OK
                    warnings.warn("Skipping exon without id: "+line)
            #elif ls[2] in ['gene', 'pseudogene']:
            elif ls[2] == 'gene':
                info['strand']=ls[6]
                info['chromosome']=chrom
                self.genes[chrom][start:end]=info
            elif all([v in info for v in ['Parent',"ID"  , 'Name']]) and info['Parent'].startswith('gene'): #those denote transcripts
                try:
                    transcripts[info["Parent"]].append((info["Name"],info["ID"]))
                except KeyError:
                    transcripts[info["Parent"]]=[(info['Name'], info['ID']) ]
            else: skipped.add(ls[2])
        print('skipped the following categories: {}'.format(skipped))
        #sort the exons
        for tid in exons.keys():
            exons[tid].sort()
        #add transcripts to genes
        for chrom in self.genes:                
            for gene in self.genes[chrom]:
                g_id=gene.data['ID']
                try:
                    t_ids=transcripts[g_id]
                except KeyError:
                    t_ids=[(gene.data['Name'],g_id)]
                for t_name, t_id in t_ids:        
                    try:
                        gene.data['transcripts'][t_name]=exons[t_id]
                    except KeyError:
                        try:
                            gene.data['transcripts']={t_name:exons[t_id]}
                        except KeyError:
                            warnings.warn('skipping {}/{} (no exons)'.format(t_name, t_id))                

    def add_splice_graphs(self, force=False):
        for chrom in self.genes:
            print('adding splice graph for chromosome '+chrom)
            for gene in tqdm(self.genes[chrom],total=len(self.genes[chrom])):
                if force or 'splice_graph' not in gene.data:
                    gene.data['splice_graph']=IntervalTree()
                    if 'transcripts' not in gene.data:
                        gene.data['splice_graph'].add(Interval(*gene[:2]), (set(), set()))
                        continue
                    for t_id,exons in gene.data['transcripts'].items():
                        for i,e in enumerate(exons):
                            pre={exons[i-1][1]} if i>0 else set()
                            suc={exons[i+1][0]} if i<len(exons)-1 else set()
                            self._add_exon_to_splice_graph(gene.data['splice_graph'], e, pre, suc)

    @staticmethod
    def _add_exon_to_splice_graph(splice_graph, exon,predecessor, successor):
        ol_exons=splice_graph.overlap(*exon)
        if not ol_exons: # no overlapping exons so far, just add to splice tree
            splice_graph.add(Interval(*exon,(predecessor, successor)))
            return
        links={exon[0]:predecessor, exon[1]:next} 
        ## the keys contain the positions, where exons start or end 
        ## the values contain the end/start positions of predecessors and successors (as a set)
        ## for exon splits, key equals one of the values
        for ol in ol_exons:
            if ol.begin in links:
                links[ol.begin].update(ol.data[0])
            else:
                links[ol.begin]=ol.data[0]
                if exon[0]<ol.begin and ol.begin <exon[1]:
                    links[ol.begin].add(ol.begin)
                else:
                    links[exon[0]].add(exon[0])
            if ol.end in links:
                links[ol.end].update(ol.data[1])
            else:
                links[ol.end]=ol.data[1]
                if exon[0]<ol.end and ol.end <=exon[1]:
                    links[ol.end].add(ol.end)
                else:
                    links[exon[1]].add(ol.end)
            splice_graph.remove(ol)                    
        links_list=sorted(links.items())
        for i,end in enumerate(links_list[1:]):
            start=links_list[i]
            splice_graph.add(Interval(start[0], end[0],
                        ({s_pos for s_pos in start[1] if s_pos <= start[0]},
                        {e_pos for e_pos in   end[1] if e_pos >=   end[0]})))

    def check_bam(self, bam_fn):
        pass

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='compare isoseq alignment to refseq')
    parser.add_argument('-i','--bam_fn', help='input bam file', type=str)
    parser.add_argument('-r','--refseq_fn', help='refseq gtf file', type=str)
    parser.add_argument('-o','--outsuffix', help='suffix for output file', default='_refseq_tab.txt')
    parser.add_argument('-a','--analyse', help='analyse table', type=str)
    args = vars(parser.parse_args())
    print(args)
    if 'bam_fn' in args and args['bam_fn'] is not None:
        args['out_fn']=os.path.splitext(args['bam_fn'])[0]+args['outsuffix']
        sg=SpliceGraph(args['refseq_fn'])
        sg.check_bam(args['bam_fn'])
    