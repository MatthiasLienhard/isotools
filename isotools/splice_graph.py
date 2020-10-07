#!/usr/bin/python3
import itertools
import warnings
from itertools import combinations
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from intervaltree import IntervalTree, Interval
from tqdm import tqdm
import logging
import scipy.stats as stats
from math import log10,pi
import isotools.transcriptome 


log=logging.getLogger(__name__)
log.setLevel(logging.INFO)
log_format=logging.Formatter('%(levelname)s: [%(asctime)s] %(name)s: %(message)s')
#log_file=logging.FileHandler('logfile.txt')
log_stream=logging.StreamHandler()
log_stream.setFormatter(log_format)
log.handlers=[] #to prevent multiple messages upon reload
log.addHandler(log_stream)

def overlap(pos1,pos2,width, height):
    if abs(pos1[0]-pos2[0])<width and abs(pos1[1]-pos2[1])<height:
        return True
    return False

class SpliceGraph():
    def __init__(self, exons, end=None, weights=None, ids=None):      
        if isinstance(exons,isotools.transcriptome.Gene):
            #end=exons.end ?uncomment to increase performance
            ids=list(exons.transcripts.keys())
            try:
                weights=np.array([tr['coverage'] for tr in exons.transcripts.values()]).swapaxes(0,1)
            except KeyError:
                weights=None
            exons=[tr['exons'] for tr in exons.transcripts.values()]
        if ids is None: #exons is a list of exons
            ids=list(range(len(exons)))
        if weights is None:
            weights=np.ones((1,len(exons)))
        self.weights=weights
        self.ids=ids
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

    def search_transcript(self, exons):
        #test if a transcript is contained in sg and returns the id(s)
        #fst special case: exons extends splice graph 
        if exons[0][1]<=self[0].start or exons[-1][0]>=self[-1].end:
            return [] 
        #snd special case: single exon transcript: return all overlapping single exon transcripts form sg
        if len(exons)==1: 
            return [self.ids[trid] for trid,(j1,j2) in enumerate(zip(self._tss,self._pas)) 
                    if self.is_same_exon(j1,j2,trid) and self[j1].start<exons[0][1] and self[j2].end > exons[0][0]]
        # all junctions must be contained and no additional 
        tr=set(range(len(self._tss)))
        j=0
        for i,e in enumerate(exons[:-1]):
            while j<len(self) and self[j].end<e[1]:#check exon (no junction allowed)
                tr -= set(trid for trid, j2 in self[j].suc.items() if self[j].end != self[j2].start)
                j+=1
            if  self[j].end!=e[1]:
                return []
            #check junction (must be present)
            tr &= set(trid for trid, j2 in self[j].suc.items() if self[j2].start == exons[i+1][0])
            j+=1
            if len(tr)==0:
                return tr

        while j<len(self):#check last exon (no junction allowed)
            tr -= set(trid for trid, j2 in self[j].suc.items() if self[j].end != self[j2].start)
            j+=1    
        return [self.ids[trid] for trid in tr]


    def is_same_exon(self, j1,j2,tr_nr):
        #test if nodes j1 and j2 belong to same exon in transcript tr_nr
        for j in range(j1,j2):
            if tr_nr not in self[j].suc or self[j].suc[tr_nr]>j+1 or self[j].end!= self[j+1].start:
                return False
        return True




    def get_alternative_splicing(self, exons,strand): 
        #returns a list of novel splicing events or splice_identical or combination
        #t='novel exon','intron retention', 'novel exonic', 'novel intronic','splice identical','combination', 'truncation', 'extention'
        tr=self.search_transcript(exons)
        if tr:
            return {'splice_identical':tr}
        is_reverse=strand=='-'
        altsplice={}
        
        j1=next(j for j,n in enumerate(self) if n.end > exons[0][0])
        #j1: index of first segment ending after exon start (i.e. first overlapping segment)        
        try:
            j2=next(j-1 for j in range(j1,len(self)) if self[j].start >= exons[0][1])
        except StopIteration:
            j2=len(self)-1
        #j2: index of last segment starting befor exon end (i.e. last overlapping segment)        

        if self[j1].pre and not any(j in self._tss for j in range(j1,j2+1)):
            tr_nr, j0=max(((trid,self._tss[trid]) for trid in self[j1].pre), key=lambda x:x[1])#j0 is the closest start node
            if not self.is_same_exon(j0,j1,tr_nr):
                end='3' if is_reverse else '5'
                altsplice.setdefault(f'{end}\' truncation',[]).append([self[j0].start, exons[0][0]]) #at start (lower position)
        for i,_ in enumerate(exons[:-1]):                               
            log.debug(f'check exon {i}:{exons[i]} between sg noded {j1}:{self[j1]} and {j2}:{self[j2]}')
            j1, exon_altsplice=self._check_exon(j1,j2,i==0,is_reverse,exons[i], exons[i+1])
            for k,v in exon_altsplice.items(): #e and the next if any
                altsplice.setdefault(k,[]).extend(v)
            if j1==len(self): #additional exons after end of splicegraph
                for remain in exons[(i+1):-1]:
                    altsplice.setdefault('novel exon',[]).append(remain)
                if i<len(exons)-1:
                    site='TSS' if is_reverse else 'PAS'
                    altsplice.setdefault(f'novel {site}',[]).append(exons[-1])
                break
            # find j2: index of last segment starting befor exon end (i.e. last overlapping  segment)
            try:
                j2=next(j-1 for j in range(j1,len(self)) if self[j].start >= exons[i+1][1])
            except StopIteration:
                j2=len(self)-1
        else:# check last exon
            j1, exon_altsplice=self._check_exon(j1,j2,False,is_reverse,exons[-1]) #todo: check is_reverse False
            for k,v in exon_altsplice.items(): #e and the next if any
                altsplice.setdefault(k,[]).extend(v)
        if self[j2].suc and not any(j in self._pas for j in range(j1,j2+1)):
            tr_nr, j3=min(((trid,self._pas[trid]) for trid in self[j2].suc), key=lambda x:x[1]) #j3 is the next end node (pas/tss on fwd/rev)
            if not self.is_same_exon(j2,j3,tr_nr):
                end='5' if is_reverse else '3'
                altsplice.setdefault(f'{end}\' truncation',[]).append([exons[-1][1],self[j3].end])
        if not altsplice:#all junctions are contained but search_transcript() did not find it
            altsplice={'novel combination':None }
        return altsplice
            
    def _check_exon(self,j1,j2,is_first,is_reverse, e, e2=None):
        #checks wether exon is supported by splice graph between nodes j1 and j2
        #j1: index of first segment ending after exon start (i.e. first overlapping segment)
        #j2: index of last segment starting befor exon end (i.e. last overlapping  segment)
        if j1>j2: #e is not contained at all  -> intronic   
            if is_first or e2==None:
                altsplice={'novel PAS' if is_first==is_reverse else 'novel TSS':[e]}
            else:
                altsplice={'novel exon':[e]}
            j2=j1
        else:
            altsplice={}
            if not is_first and self[j1][0]!=e[0]:
                pos="intronic" if self[j1][0]>e[0] else "exonic"
                kind='donor' if is_reverse else 'acceptor'
                dist=self[j1][0]-e[0] if (e[0]-self[j1][0]<self[j1][1]-e[0]) else self[j1][1]-e[0] #the smaller of the two distances to next junction
                altsplice[f'novel {pos} splice {kind}']=[(e[0],dist)] 
            if e2 is not None and self[j2][1]!=e[1]:
                pos="intronic" if self[j2][1]<e[1] else "exonic"
                kind='acceptor' if is_reverse else 'donor'
                dist= self[j2][0]-e[1] if (e[1]-self[j2][0]<self[j2][1]-e[1]) else self[j2][1]-e[1] 
                altsplice.setdefault(f'novel {pos} splice {kind}',[]).append((e[1],dist))
            #find retained introns
            for j in range(j1,j2):
                if self[j].suc:
                    gap,j_suc=min(((self[j_suc][0]-self[j][1],j_suc) for j_suc in set(self[j].suc.values()) ), key=lambda x: x[0])
                    #todo: (potentially pedantic note) check transcript supporting this intron should start befor e[1] - 
                    # otherwise alternative tss within the retained intron might give unwanted results
                    if gap>0 and self[j_suc][0] < e[1]:                         
                        altsplice.setdefault('retained intron',[]).append([self[j][1],self[j_suc][0]])   
            #check for overlapping retained introns, take the minimum length
            if 'retained intron' in altsplice and len(altsplice['retained intron'])>1:
                ret_introns=altsplice['retained intron']
                rm=[]
                for idx,intron in enumerate(ret_introns):
                    if idx in rm:
                        continue
                    ol=[idx]
                    idx2=idx+1
                    while idx2<len(ret_introns) and ret_introns[idx2][0]<= intron[1]:
                        if idx2 not in rm:
                            ol.append(idx2)
                        idx2+=1
                    if len(ol)>1:  #several overlapping retained introns: keep only the smallest
                        keep=min(ol, key=lambda x: ret_introns[x][1]-ret_introns[x][0])
                        rm.extend([x for x in ol if x != keep])
                altsplice['retained intron']=[intron for idx,intron in enumerate(ret_introns) if idx not in rm]
                #             
            j1=j2 #j1 becomes index of last segment starting befor e ends
        if e2 is not None: #check presence of junction
            introns=[]
            while j2 < len(self) and self[j2][1]<=e2[0]:
                if j2 in self[j2-1].suc.values() and self[j2-1][1]<self[j2][0]:
                    introns.append([self[j2-1][1],self[j2][0]])
                j2+=1    #j2 now is index of first segment ending after e2 starts (first overlapping)
            if j2 in self[j2-1].suc.values() and self[j2-1][1]<self[j2][0]:
                introns.append([self[j2-1][1],self[j2][0]])    
            if j2<len(self) and self[j1][1]==e[1] and self[j2][0]==e2[0] and j2 not in self[j1].suc.values():#junciton not present...
                if len(introns)>1: #strictly we would need to check the path from j1 to j2
                    for i in range(len(introns)-1):
                        altsplice.setdefault('exon skipping',[]).append([introns[i][1], introns[i+1][0]])
                else:
                    altsplice.setdefault('novel junction',[]).append([e[1],e2[0]]) #for example mutually exclusive exons spliced togeter
        log.debug(f'check exon {e} resulted in {altsplice}')
        return j2, altsplice
            
    def fuzzy_junction(self,exons,size):
        #for each intron from "exons", look for introns in the splice graph shifted by less than "size".
        # these shifts may be produced by ambigious alignments.
        #return a dict with the intron number as key and the shift as value 
        #assuming size is smaller than introns
        fuzzy={}
        if size < 1: #no need to check
            return fuzzy
        j1=0
        for i,e1 in enumerate(exons[:-1]):
            e2=exons[i+1]
            try:
                j1=next(j for j in range(j1,len(self)) if self[j].end+size >=e1[1] )
                #j1: first node intersecting size range of e1 end
            except StopIteration:
                break
            
            shift=None
            while j1 < len(self) and self[j1].end - e1[1] <= size: # in case there are several nodes starting in the range around e1
                shift_e1=self[j1].end - e1[1]
                #print(f'{i} {e1[1]}-{e2[0]} {shift_e1}')
                if shift_e1==0:
                    break
                if any( self[j2].start - e2[0] == shift_e1 for j2 in set(self[j1].suc.values()) ):
                    shift=shift_e1
                j1+=1
            else:
                if shift is not None:
                    fuzzy[i]=shift
        return fuzzy

    def get_intersects(self, exons):
        #returns the splice junction and exonic base intersects
        intersect=[0,0]
        i=j=0
        while True:            
            if self[j][0] == exons[i][0] and any(self[k][1] < self[j][0] for k in self[j].pre.values()):
                intersect[0]+=1 # same position and actual splice junction(not just tss or pas and internal junction)
            if self[j][1] == exons[i][1] and any(self[k][0] > self[j][1] for k in self[j].suc.values()):
                intersect[0]+=1
            if self[j][1] > exons[i][0] and exons[i][1] > self[j][0]: #overlap
                intersect[1]+=min(self[j][1], exons[i][1])-max(self[j][0], exons[i][0])
            if exons[i][1]<self[j][1]:
                i+=1
            else:
                j+=1
            if i==len(exons) or j==len(self):
                return(intersect)


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
                            if pvalue<pval_th:
                                #print(f'found {i}:{starts[i]} vs {j}:{ends[j]}\n{ma}\n:p={pvalue}\n-------')
                                alt_from=['tss' if k == -1 else self[k].end for k,_ in [starts[i][idx] for idx in ixpair]]
                                alt_to=['pas' if k == -1 else self[k].start for k,_ in [ends[j][idy] for idy in iypair]]
                                found.append((pvalue,subma, alt_from,self[i].start,self[j].end,alt_to))
        return found

    def get_splice_coverage(self):
        '''search for alternative paths in the splice graphs ("loops"), 
        e.g. combinations of nodes A and B with more than one path from A to B'''
        for i, nA in enumerate(self):
            junctions=sorted(list({j for j in nA.suc.values()})) # target nodes for junctions from node A ordered by intron size
            if len(junctions)<2:
                continue #no alternative
            j_sets={}#transcripts supporting the different  junctions
            for tr,node_id in nA.suc.items(): 
                j_sets.setdefault(node_id,[]).append(tr)
            log.debug(f'checking node {i}: {nA} ({junctions})')
            for idx,joi in enumerate(junctions[1:]): # start from second, as first does not have an alternative
                nB=self[joi]
                #alternative=[tr for j in junctions[:idx+1] for tr in j_sets[j] if self._pas[tr] > joi]
                alternative={j:[tr for tr in j_sets[j] if self._pas[tr] > joi] for j in junctions[:idx+1] }
                log.debug(f'loop for {joi} from {nA.end} to {nB.start}: joi={j_sets[joi]}, alternatives={alternative}')
                if any(a for a in alternative.values()):   
                    #get type of alternative splicing                 
                    alt_type=set()
                    for j,j_set in alternative.items():
                        if self[j].start==nA.end:
                            if self[j].end==nB.start:
                                alt_type.add('IR')#intron retention
                            else:
                                alt_type.add('AS1')#alternative splice site (5' on +, 3' on -)
                        elif j+1==joi:
                            if self[j].end==nB.start:
                                alt_type.add('AS2')#alternative splice site (5' on -, 3' on +)
                            else:
                                alt_type.add('ES')#exon skiping
                        else:
                            #more than one node
                            alt_type.add('UK')
                        if not all(j in nB.pre for j in j_set):
                            alt_type.add('ME')#mutually exclusive exon/sequence
                    weight=self.weights[:,j_sets[joi]].sum(1)                            
                    alt_weight=self.weights[:,[a for l in alternative.values() for a in l]].sum(1)
                    yield weight, weight+alt_weight,nA.end,nB.start,alt_type
                    
                


    def get_splice_coverage_old(self):
        for i,first in enumerate(self):
            for j,second in ((idx,self[idx]) for idx in set(first.suc.values())):
                if first.end== second.start: #only consider splice junctions (this excludes alternative tss and pas... todo)
                    continue
                idx=[tr for tr,idx in first.suc.items() if idx==j] #transcripts supporting this junction
                relevant=[idx for idx,(tss,pas) in enumerate(zip(self._tss, self._pas)) if tss<=i and pas>=j] #transcripts spanning this position
                if len(idx)==len(relevant): #no additional spanning transcripts
                    continue
                total_weight=self.weights[:,relevant].sum(1)
                weight=self.weights[:,idx].sum(1)
                #how to define the type (e.g. exon skipping/intron retention/alternative 5'/3')
                yield weight,total_weight,first.end,second.start




        
    def __getitem__(self, key):
        return self._graph[key]
    
    def __len__(self):
        return len(self._graph)

    #def pre(self,key):
    #    return (self._graph[pre_key] for pre_key in self._graph[key][2])
    #    
    #def suc(self,key):
    #    if self._graph[key][3]:
    #        return (self._graph[suc_key] for suc_key in self._graph[key][3])
    #    else:
    #        return self._graph[key][3]
    
    
        

class SpliceGraphNode(tuple):
    def __new__(cls,start, end, pre=None, suc=None):
        if pre is None:
            pre=dict()
        if suc is None:
            suc=dict()
        return super(SpliceGraphNode,cls).__new__(cls,(start,end,pre, suc))
    def __getnewargs__(self):
        return tuple(self)
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
