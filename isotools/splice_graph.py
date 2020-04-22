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
import matplotlib.patches as patches

from intervaltree import IntervalTree, Interval
from tqdm import tqdm
import logging
import scipy.stats as stats
from math import log10,pi
import isotools.transcriptome 
from scipy.stats import binom, chi2

def overlap(pos1,pos2,width, height):
    if abs(pos1[0]-pos2[0])<width and abs(pos1[1]-pos2[1])<height:
        return True
    return False

class SpliceGraph():
    def __init__(self, exons, end=None, weights=None):      
        if isinstance(exons,isotools.transcriptome.Gene):
            #end=exons.end ?uncomment to increase performance
            self.ids=list(exons.transcripts.keys())
            weights=np.array([tr['coverage'] for tr in exons.transcripts.values()]).swapaxes(0,1)
            exons=[tr['exons'] for tr in exons.transcripts.values()]
        else: #exons is a list of exons
            self.ids=list(range(exons))
        if weights is None:
            weights=np.ones((1,len(exons)))
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
        start_idx={node.start:i for i,node in enumerate(self._graph)}
        end_idx={node.end:i for i,node in enumerate(self._graph)}
        
        self._tss=[start_idx[e[0][0]] for e in exons]
        self._pas=[end_idx[e[-1][1]] for e in exons]

        for i,tr in enumerate(exons):
            for j,e in enumerate(tr):
                if start_idx[e[0]]<end_idx[e[1]]: # psydojuctions within exon
                    self._graph[start_idx[e[0]]].suc[i]=start_idx[e[0]]+1
                    self._graph[end_idx[e[1]]].pre[i]=end_idx[e[1]]-1
                    for node_idx in range(start_idx[e[0]]+1,end_idx[e[1]]):
                        self._graph[node_idx].suc[i]=node_idx+1
                        self._graph[node_idx].pre[i]=node_idx-1
                if j<len(tr)-1:# real junctions
                    e2=tr[j+1]
                    self._graph[end_idx[e[1]]].suc[i]=start_idx[e2[0]]
                    self._graph[start_idx[e2[0]]].pre[i]=end_idx[e[1]]
                    
    def restore(self, i):    #mainly for testing    
        idx=self._tss[i]
        exons=[[self._graph[idx].start, self._graph[idx].end]]
        while True:
            idx=self._graph[idx].suc[i]
            if self._graph[idx].start==exons[-1][1]:#extend
                exons[-1][1]=self._graph[idx].end
            else:
                exons.append([self._graph[idx].start, self._graph[idx].end])
            #print(exons)
            if idx == self._pas[i]:
                break
        return exons

    def restore_reverse(self, i):       #mainly for testing     
        idx=self._pas[i]
        exons=[[self._graph[idx].start, self._graph[idx].end]]
        while True:
            idx=self._graph[idx].pre[i]
            if self._graph[idx].end==exons[-1][0]:#extend
                exons[-1][0]=self._graph[idx].start
            else:
                exons.append([self._graph[idx].start, self._graph[idx].end])
            #print(exons)
            if idx == self._tss[i]:
                break
        exons.reverse()
        return exons

    def find_ts_candidates(self):        
        for i, gnode in enumerate(self._graph[:-1]):
            if self._graph[i+1].start==gnode.end: #jump candidates: introns that start within an exon
                jumps={idx:n for idx,n in gnode.suc.items() if n>i+1 and self._graph[n].start==self._graph[n-1].end}
                #find jumps (n>i+1) and check wether they end within an exon begin(jumptarget)==end(node before)
                jump_weight={}
                for idx, target in jumps.items():
                    jump_weight.setdefault(target, [0,[]])
                    jump_weight[target][0]+=self.weights[:,idx].sum(0)
                    jump_weight[target][1].append(idx)

                for target, (w,idx) in jump_weight.items():
                    long_idx=set(idx for idx,n in gnode.suc.items() if n==i+1) & set(idx for idx,n in self[target].pre.items() if n==target-1)
                    try:
                        longer_weight=self.weights[:,list(long_idx)].sum()#sum(self.weights[i] for i in long_idx)
                    except IndexError:
                        print(long_idx)
                        raise
                    yield gnode.end, self[target].start,w, longer_weight, [self.ids[i] for i in idx]
        
    def splice_dependence(self, pval_th=0.05, min_sum=10):
        #starts[i]: list with paths that join at node i
        #end[i]: list with paths that split at node i
        # the lists contain tupels (i, [trids]) where i is the neighbour exon number (-1 at tss/pas) and a list of transcript numbers
        n=len(self)
        starts=[list() for _ in range(n)]
        ends=[list() for _ in range(n)]
        for i, gnode in enumerate(self):
            if i in self._tss: 
                starts[i].append((-1, [j for j,v in enumerate(self._tss) if v==i]))
            if i in self._pas: 
                ends[i].append((-1, [j for j,v in enumerate(self._pas) if v==i]))
            starts[i].extend([(j, [k for k in gnode.pre if gnode.pre[k]==j]) for j in set(gnode.pre.values())])
            ends[i].extend([(j, [k for k in gnode.suc if gnode.suc[k]==j]) for j in set(gnode.suc.values())])
        found=list()
        #test independence of starts[i] and ends[j] for i,j where i>=j
        for i in range(len(starts)):
            for j in range(i,len(ends)):
                if len(starts[i])>1 and len(ends[j])>1:#there are options at both i and j   
                    ma=np.zeros((len(starts[i]),len(ends[j]))) #the contingency table
                    for ix, st in enumerate(starts[i]):
                        for iy, en in enumerate(ends[j]):
                            ma[ix][iy]=self.weights[:,[k for k in st[1] if k in en[1]]].sum()  #sum the weight for the test
                    
                    #fisher test requires 2x2 tables, and filter by row and colum sum
                    ix=[i for i,s in enumerate(ma.sum(axis=1)) if s > min_sum]
                    iy=[i for i,s in enumerate(ma.sum(axis=0)) if s > min_sum]
                    for ixpair in combinations(ix, 2):
                        for iypair in combinations(iy, 2):     
                            subma=  ma[np.ix_(ixpair, iypair)]
                            if any(subma.sum(0)<min_sum) or any(subma.sum(1)<min_sum):
                                continue
                            oddsratio, pvalue = stats.fisher_exact(subma)
                    #fisher test requires 2x2 tables, and filter by row and colum sum
                    ix=[i for i,s in enumerate(ma.sum(axis=1)) if s > min_sum]
                    iy=[i for i,s in enumerate(ma.sum(axis=0)) if s > min_sum]
                    for ixpair in combinations(ix, 2):
                        for iypair in combinations(iy, 2):     
                            subma=  ma[np.ix_(ixpair, iypair)]
                            if any(subma.sum(0)<min_sum) or any(subma.sum(1)<min_sum):
                                continue
                            oddsratio, pvalue = stats.fisher_exact(subma)
                            if pvalue<pval_th:
                                #print(f'found {i}:{starts[i]} vs {j}:{ends[j]}\n{ma}\n:p={pvalue}\n-------')
                                alt_from=['tss' if k == -1 else self[k].end for k,_ in [starts[i][idx] for idx in ixpair]]
                                alt_to=['pas' if k == -1 else self[k].start for k,_ in [ends[j][idy] for idy in iypair]]
                                found.append((pvalue,subma, alt_from,self[i].start,self[j].end,alt_to))
        return found

    def altsplice_test(self, groups,min_junction_cov=10):      #todo: considier min_junction_cov  
        assert len(groups)==2
        assert all(isinstance(g,list) for g in groups)
        p=dict()
        #total_weight=[self.weights[grp,:].sum(1) for grp in groups]
        for i,(es,ee, pre, suc) in enumerate(self):
            for target in set(suc.values()):
                if self[target][0]== ee: #only consider splice junctions (this excludes alternative tss and pas... todo)
                    continue
                relevant=[i for i,(tss,pas) in enumerate(zip(self._tss, self._pas)) if tss<=i and pas>=target]
                total_weight=[self.weights[np.ix_(grp,relevant)].sum(1) for grp in groups]
                weight=[self.weights[np.ix_(grp,[tr for tr in suc if suc[tr]==target])].sum(1) for grp in groups]
                #proportion estimates
                phat_0=sum(w+1 for wg in weight for w in wg)/sum(tw+2 for twg in total_weight for tw in twg) #null hypothesis (just one proportion)
                phat_1=[sum(w+1 for w in wg)/sum(tw+2 for tw in twg) for wg,twg in zip(weight, total_weight)] #alternative: one proportion per group
                # calculate the log likelihoods
                l0 = binom.logpmf(weight, total_weight, phat_0).sum()
                l1 = binom.logpmf(weight,total_weight,[[p]*len(groups[i]) for i,p in enumerate(phat_1)]).sum()
                # calculate the pvalue (sf=1-csf(), 1df)
                p[(ee,self[target][0])]=[chi2.sf(2*(l1-l0),1)]+phat_1+[tw.sum() for tw in total_weight]
        return p

    
                
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
        start_ref={e.start:i for i,e in enumerate(self)}
        end_ref={e.end:i for i,e in enumerate(self)}
        
        merge_to=[[i] for i,_ in enumerate(self)]
        for i,node in enumerate(self):
            if node.end-node.start < fuzzy_junction: # short psydoexon
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
                starts=[self[i].start for i in merge_set]
                ends=[self[i].end for i in merge_set]
                pos=(min(starts), max(ends)) #TODO: CURRENTLY WE TAKE THE OUTER, BUT JUNCTION with best support WOULD BE BETTER?
                new_start.update({self[i].start:pos[0] for i in merge_set})
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
    def __new__(cls,start, end, pre=None, suc=None):
        if pre is None:
            pre=dict()
        if suc is None:
            suc=dict()
        return super(SpliceGraphNode,cls).__new__(cls,(start,end,pre, suc))
    
    @property
    def start(self):
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
