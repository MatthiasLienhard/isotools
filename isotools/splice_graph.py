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
import logging

import isotools.transcriptome



class SpliceGraph():
    def __init__(self, exons, end=None, weights=None):      
        if isinstance(exons,isotools.transcriptome.Gene):
            weights=[tr['nZMW'] for tr in exons.transcripts]
            exons=[tr['exons'] for tr in exons.transcripts]
            
        self.weights=weights
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
        
        self._tss=[begin_idx[e[0][0]] for e in exons]
        self._pas=[end_idx[e[-1][1]] for e in exons]

        for i,tr in enumerate(exons):
            for j,e in enumerate(tr):
                if begin_idx[e[0]]<end_idx[e[1]]: # psydojuctions within exon
                    self._graph[begin_idx[e[0]]].suc[i]=begin_idx[e[0]]+1
                    self._graph[end_idx[e[1]]].pre[i]=end_idx[e[1]]-1
                    for node_idx in range(begin_idx[e[0]]+1,end_idx[e[1]]):
                        self._graph[node_idx].suc[i]=node_idx+1
                        self._graph[node_idx].pre[i]=node_idx-1
                if j<len(tr)-1:# real junctions
                    e2=tr[j+1]
                    self._graph[end_idx[e[1]]].suc[i]=begin_idx[e2[0]]
                    self._graph[begin_idx[e2[0]]].pre[i]=end_idx[e[1]]
                    
    def restore(self, i):    #mainly for testing    
        idx=self._tss[i]
        exons=[[self._graph[idx].begin, self._graph[idx].end]]
        while True:
            idx=self._graph[idx].suc[i]
            if self._graph[idx].begin==exons[-1][1]:#extend
                exons[-1][1]=self._graph[idx].end
            else:
                exons.append([self._graph[idx].begin, self._graph[idx].end])
            #print(exons)
            if idx == self._pas[i]:
                break
        return exons

    def restore_reverse(self, i):       #mainly for testing     
        idx=self._pas[i]
        exons=[[self._graph[idx].begin, self._graph[idx].end]]
        while True:
            idx=self._graph[idx].pre[i]
            if self._graph[idx].end==exons[-1][0]:#extend
                exons[-1][0]=self._graph[idx].begin
            else:
                exons.append([self._graph[idx].begin, self._graph[idx].end])
            #print(exons)
            if idx == self._tss[i]:
                break
        exons.reverse()
        return exons

    def find_ts_candidates(self):
        
        for i, gnode in enumerate(self._graph[:-1]):
            if self._graph[i+1].begin==gnode.end: #jump candidates: introns that start within an exon
                jumps={idx:n for idx,n in gnode.suc.items() if n>i+1 and self._graph[n].begin==self._graph[n-1].end}
                #find jumps (n>i+1) and check wether they end within an exon begin(jumptarget)==end(node before)
                jump_weight={}
                for idx, target in jumps.items():
                    jump_weight.setdefault(target, [0,[]])
                    jump_weight[target][0]+=self.weights[idx]
                    jump_weight[target][1].append(idx)

                for target, (w,idx) in jump_weight.items():
                    long_idx=set(idx for idx,n in gnode.suc.items() if n==i+1) & set(idx for idx,n in self[target].pre.items() if n==target-1)
                    longer_weight=sum(self.weights[i] for i in long_idx)
                    yield gnode.end, self[target].begin,w, longer_weight, idx
                    


    '''
    def find_5truncations(self, is_reverse=False, dist_3=10):
        truncations=list()
        next_idx=2 if is_reverse else 3
        grps=dict()
        start_idx=self._tss if is_reverse else self._pas
        for i,v in enumerate(start_idx):
            grps.setdefault(v, set()).add(i)
        for idx, grp in grps.items():
            pass
            #
        return truncations
    ''' 
        
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
    
    
    
    '''
    def unify_junctions(self, fuzzy_junction=5, reference=None):
        start_ref={e.begin:i for i,e in enumerate(self)}
        end_ref={e.end:i for i,e in enumerate(self)}
        
        merge_to=[[i] for i,_ in enumerate(self)]
        for i,node in enumerate(self):
            if node.end-node.begin < fuzzy_junction: # short psydoexon
                pre=node.pre #todo: currently merge with predecessor has prioritiy, but considering junction support would be better
                if len(pre)==1 and self[pre[0]].end == node.start:#merge with predesessor                
                    merge_to[node.pre]+=merge_to[i]
                    merge_to[i]=merge_to[node.pre]
                else:
                    suc=node.suc
                    if len(suc)==1 and self[suc[0]].start ==node.end:#merge with successor
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
                starts=[self[i].begin for i in merge_set]
                ends=[self[i].end for i in merge_set]
                pos=(min(starts), max(ends)) #TODO: CURRENTLY WE TAKE THE OUTER, BUT JUNCTION with best support WOULD BE BETTER?
                new_start.update({self[i].begin:pos[0] for i in merge_set})
                new_end.update({self[i].end:pos[1] for i in merge_set})
                pre={p for i in merge_set for p in self[i].pre if p not in starts}
                suc={s for i in merge_set for s in self[i].suc if s not in ends}
                graph.append(SpliceGraphNode(*pos,pre, suc))
                if any(idx in self._tss for idx in merge_set):
                    tss.append(i)
                if any(idx in self._pas for idx in merge_set):
                    pas.append(i)            
                i+=1
        self._graph =graph
        self._tss=tss
        self._pas=pas
        return new_start, new_end
    '''
class SpliceGraphNode(tuple):
    def __new__(cls,begin, end, pre=None, suc=None):
        if pre is None:
            pre=dict()
        if suc is None:
            suc=dict()
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
