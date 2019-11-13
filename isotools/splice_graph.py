#!/usr/bin/python3
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

import isotools.utils
import isotools.collapse
                          


                
def get_splice_graph(gene):
    return(SpliceGraph((tr['exons'] for tr in gene.data['transcripts'].values()), gene.end))


def unify_junctions_of_gene(gene, fuzzy_junction=5):
    if 'splice_graph' in gene.data:
        splice_graph=gene.data['splice_graph']
    else:
        splice_graph=isotools.splice_graph.get_splice_graph(gene)
        gene.data['splice_graph'] =splice_graph

    start_ref={e.begin:i for i,e in enumerate(splice_graph)}
    end_ref={e.end:i for i,e in enumerate(splice_graph)}
    merge_to=[[i] for i,_ in enumerate(splice_graph)]
    for i,node in enumerate(splice_graph):
        if node.end-node.begin < fuzzy_junction: # short psydoexon
            pre=node.pre #todo: currently merge with predecessor has prioritiy, but considering junction support would be better
            if len(pre)==1 and splice_graph[pre[0]].end == node.start:#merge with predesessor                
                merge_to[node.pre]+=merge_to[i]
                merge_to[i]=merge_to[node.pre]
            else:
                suc=node.suc
                if len(suc)==1 and splice_graph[suc[0]].start ==node.end:#merge with successor
                    merge_to[node.suc]+=merge_to[i]
                    merge_to[i]=merge_to[node.suc]
    seen=set()
    graph=list()
    new_start={}
    new_end={}
    i=0
    tss=[]
    pas=[]
    for merge_set in merge_to:
        if id(merge_set) not in seen:
            seen.add(id(merge_set))
            starts=[splice_graph[i].begin for i in merge_set]
            ends=[splice_graph[i].end for i in merge_set]
            pos=(min(starts), max(ends)) #TODO: CURRENTLY WE TAKE THE OUTER, BUT JUNCTION with best support WOULD BE BETTER?
            new_start.update({splice_graph[i].begin:pos[0] for i in merge_set})
            new_end.update({splice_graph[i].end:pos[1] for i in merge_set})
            pre={p for i in merge_set for p in splice_graph[i].pre if p not in starts}
            suc={s for i in merge_set for s in splice_graph[i].suc if s not in ends}
            graph.append(SpliceGraphNode(*pos,pre, suc))
            if any(idx in splice_graph._tss for idx in merge_set):
                tss.append(i)
            if any(idx in splice_graph._pas for idx in merge_set):
                pas.append(i)            
            i+=1
    splice_graph._graph =graph
    splice_graph._tss=tss
    splice_graph._pas=pas
    for n,tr in gene.data['transcripts'].items():
        try:
            tr['exons']=[(new_start[e[0]], new_end[e[1]]) for e in tr['exons']]
        except KeyError:
            print('no new values for transcript {}: {}\nnew starts:{}\nnew ends:{}'.format(n,tr['exons'],new_start, new_end))
            raise
    return gene




class SpliceGraph():
    def __init__(self, exons, end=None):        
        if end is None:
            end=max((e[-1][1] for e in exons))
        open_exons=dict()
        for tr in exons:
            for e in tr:
                open_exons[e[0]]=open_exons.setdefault(e[0],0)+1
                open_exons[e[1]]=open_exons.setdefault(e[1],0)-1
        #sort by value
        boundaries=sorted(open_exons)
        open_count=0
        self._graph=list()
        for i,start in enumerate(boundaries[:-1]):
            open_count+=open_exons[start]
            if open_count>0:
                self._graph.append(SpliceGraphNode(start, boundaries[i+1]))

        #get the links
        begin_idx={node.begin:i for i,node in enumerate(self._graph)}
        end_idx={node.end:i for i,node in enumerate(self._graph)}

        for tr in exons:
            i=0
            for node in self._graph:
                while node.begin>tr[i][1]:#next exon
                    i+=1
                    if i==len(tr): break #no more exons - break both loops
                if i==len(tr): break 
                if node.end <= tr[i][0]:
                    continue #next node
                if tr[i][0]==node.begin and i>0: #start match
                    node.pre.add(end_idx[tr[i-1][1]]) #predesessor (real splice junction)
                elif tr[i][0]<node.begin: #start overlap
                    node.pre.add(end_idx[node.begin]) #predesessor (psydojunktion within exon)
                if tr[i][1]==node.end and i<len(tr)-1: #end match
                    node.suc.add(begin_idx[tr[i+1][0]]) #succsessor
                elif tr[i][1]>node.end: #end overlap
                    node.suc.add(begin_idx[node.end]) #successor

        self._tss=[begin_idx[e[0][0]] for e in exons]
        self._pas=[end_idx[e[-1][1]] for e in exons]
        
    def __getitem__(self, key):
        return self._graph[key]
    
    def __len__(self):
        return len(self._graph)

    def pre(self,key):
        return (self._graph[pre_key] for pre_key in self._graph[key][2])
        
    def suc(self,key):
        if self._graph[key][3]:
            return (self._graph[suc_key] for suc_key in self._graph[key][3])
        else:
            return self._graph[key][3]
    
class SpliceGraphNode(tuple):
    def __new__(cls,begin, end, pre=None, suc=None):
        if pre is None:
            pre=set()
        if suc is None:
            suc=set()
        return super(SpliceGraphNode,cls).__new__(cls,(begin,end,pre, suc))
    
    @property
    def begin(self):
        return self.__getitem__(0)
    @property
    def end(self):
        return self.__getitem__(1)
    @property
    def pre(self):
        return self.__getitem__(2)
    @property
    def suc(self):
        return self.__getitem__(3)

'''
for tr in exons:
            i=0
            for j,e in enumerate(chain):
                while chain[j][0]>tr[i][1]:#chain is ahead
                    i+=1
                    if i==len(tr): break #2
                if i==len(tr): break 
                if tr[i][0]==e[0] and i>0: #start match
                    nodes[j][0].add(tr[i-1][1]) #predesessor (junction)
                elif tr[i][0]<e[0]: #start overlap
                    nodes[j][0].add(e[0]) #predesessor (within exon)
                if tr[i][1]==e[1] and i<len(tr)-1: #end match
                    nodes[j][1].add(tr[i+1][0]) #succsessor
                elif tr[i][1]>e[1]: #end overlap
                    nodes[j][1].add(e[1]) #successor
        self._graph=[SpliceGraphNode(*pos, *link) for pos,link in zip(chain,nodes)]
'''

if __name__ == '__main__':
    
    arg_parser = argparse.ArgumentParser(
        description='compare isoseq alignment to refseq')
    arg_parser.add_argument('-i', '--bam_fn',    help='input bam file', type=str, default='/project/42/pacbio/hecatos_isoseq/05-gmap/all_merged_hq_isoforms.bam')
    arg_parser.add_argument('-r', '--refseq_fn', help='refseq gtf file', type=str)
    arg_parser.add_argument('-o', '--outsuffix', help='suffix for output file', type=str, default='_refseq_tab.txt')
    arg_parser.add_argument('-a', '--analyse',   help='analyse table', type=str)
    args=arg_parser.parse_args()
    #isoseq=import_bam_transcripts(args.bam_fn, n_lines=1000)
    #collapse_transcripts(isoseq)
