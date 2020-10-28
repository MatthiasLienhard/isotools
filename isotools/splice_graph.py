from itertools import combinations
import pandas as pd
import numpy as np
from tqdm import tqdm
import logging
import scipy.stats as stats
import logging

class SpliceGraph():
    def __init__(self, transcripts, samples=None, weights=None):      
        if weights is None:
            self.weights=np.ones(len(transcripts)) #todo: is this necesary?
        else:
            self.weights=weights
        open_exons=dict()
        for tr in transcripts:
            for e in tr:
                open_exons[e[0]]=open_exons.get(e[0],0)+1
                open_exons[e[1]]=open_exons.get(e[1],0)-1
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
        self._tss=[start_idx[e[0][0]] for e in transcripts] #todo: this is missleading: on - strand this would not be the tss
        self._pas=[end_idx[e[-1][1]] for e in transcripts]
            

        for i,tr in enumerate(transcripts):
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
            if idx == self._pas[i]:
                break
            idx=self._graph[idx].suc[i]
            if self._graph[idx].start==exons[-1][1]:#extend
                exons[-1][1]=self._graph[idx].end
            else:
                exons.append([self._graph[idx].start, self._graph[idx].end])
            #print(exons)
            
        return exons

    def restore_reverse(self, i):       #mainly for testing     
        idx=self._pas[i]
        exons=[[self._graph[idx].start, self._graph[idx].end]]
        while True:
            if idx == self._tss[i]:
                break
            idx=self._graph[idx].pre[i]
            if self._graph[idx].end==exons[-1][0]:#extend
                exons[-1][0]=self._graph[idx].start
            else:
                exons.append([self._graph[idx].start, self._graph[idx].end])
            #print(exons)
            
        exons.reverse()
        return exons

    def search_transcript(self, exons):
        #test if a transcript is contained in sg and returns the id(s)
        #fst special case: exons extends splice graph 
        if exons[0][1]<=self[0].start or exons[-1][0]>=self[-1].end:
            return [] 
        #snd special case: single exon transcript: return all overlapping single exon transcripts form sg
        if len(exons)==1: 
            return [trid for trid,(j1,j2) in enumerate(zip(self._tss,self._pas)) 
                    if self.is_same_exon(j1,j2,trid) and self[j1].start<=exons[0][1] and self[j2].end >= exons[0][0]]
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
        return [trid for trid in tr]


    def is_same_exon(self, j1,j2,tr_nr):
        #test if nodes j1 and j2 belong to same exon in transcript tr_nr
        for j in range(j1,j2):
            if tr_nr not in self[j].suc or self[j].suc[tr_nr]>j+1 or self[j].end!= self[j+1].start:
                return False
        return True


    def find_truncations(self):
        truncated=set()
        contains={}
        starts=[set() for _ in range(len(self))]
        for trid,idx in enumerate(self._tss):
            starts[idx].add(trid)
        for trid,idx in enumerate(self._tss):
            if trid in truncated:
                continue
            start_in=set() #transcripts that are contained until idx
            contains[trid]=set()
            while True:
                start_in.update(starts[idx])
                contains[trid].update({c for c in start_in if self._pas[c]==idx}) # add fully contained
                if idx == self._pas[trid]:
                    contains[trid].remove(trid) #remove self
                    truncated.update(contains[trid]) #those are not checkted 
                    break
                suc=self._graph[idx].suc
                start_in={c for c in start_in if c in suc and suc[c] ==suc[trid] } #remove transcripts that split at idx
                idx=suc[trid]
        truncations={}
        for big,smallL in contains.items():
            if big not in truncated:
                for trid in smallL:
                    truncations.setdefault(trid,[]).append(big)
        return truncations


    def get_alternative_splicing(self, exons,strand): 
        'compares exons to splice graph and returns list of novel splicing events'
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
            logging.debug(f'check exon {i}:{exons[i]} between sg noded {j1}:{self[j1]} and {j2}:{self[j2]}')
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
        #checks whether exon is supported by splice graph between nodes j1 and j2
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
        logging.debug(f'check exon {e} resulted in {altsplice}')
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
                    yield gnode.end, self[target].start,w, longer_weight,  idx
        
    def splice_dependence(self, sidx,min_cov):
        tr_cov=self.weights[sidx,:].sum(0)
        if tr_cov.sum()<2* min_cov:
            return #no chance of getting above the coverage
        jumps={}
        for i,n in enumerate(self):
            for trid,suc_i in n.suc.items():
                if suc_i>i+1:
                    jumps.setdefault((i,suc_i),set()).add(trid)
        nonjumps={}
        for j in jumps:
            if tr_cov[list(jumps[j])].sum()<min_cov:
                continue
            jstart={trid for trid,suc_i in self[j[0]].suc.items() if suc_i<j[1]}
            jend={trid for trid,pre_i in self[j[1]].pre.items() if pre_i>j[0]}
            nojump_ids=jstart.intersection(jend)
            if tr_cov[list(nojump_ids)].sum()>=min_cov:
                nonjumps[j]=jstart.intersection(jend)
        for j1 in nonjumps:
            for j2 in (j for j in nonjumps if j[0]>j1[1]):
                ma=np.array([
                    tr_cov[list(jumps[j1].intersection(jumps[j2]))].sum(),
                    tr_cov[list(jumps[j1].intersection(nonjumps[j2]))].sum(),
                    tr_cov[list(nonjumps[j1].intersection(jumps[j2]))].sum(),
                    tr_cov[list(nonjumps[j1].intersection(nonjumps[j2]))].sum()]).reshape((2,2))
                if any(ma.sum(0)<min_cov) or any(ma.sum(1)<min_cov):
                    continue
                #oddsratio, pvalue = fisher_exact(ma)                
                #log.debug('j1={}, j2={}, ma={}, or={}, p={}'.format(j1,j2,list(ma),oddsratio, pvalue))
                yield ma,(self[j1[0]].end, self[j1[1]].start), (self[j2[0]].end, self[j2[1]].start)


    def get_splice_coverage(self):
        '''search for alternative paths in the splice graphs ("loops"), 
        e.g. combinations of nodes A and B with more than one path from A to B.
        returns the coverage of most direkt path and alternative paths, positions and type of alternative'''
        alt_types=['ES','ASL', 'ASR','IR','ME'] # alternative types: intron retention, alternative splice site at left and right, exon skipping, mutually exclusive
        inB_sets=[(set(),set())] #list of spliced and unspliced transcripts joining in B
        for i, nB in enumerate(self[1:]):
            inB_sets.append((set(),set()))
            unspliced=self[i].end==nB.start
            for tr,node_id in nB.pre.items():
                inB_sets[i+1][unspliced and node_id==i].add(tr)
        for i, nA in enumerate(self):
            junctions=sorted(list(set(nA.suc.values()))) # target nodes for junctions from node A ordered by intron size
            if len(junctions)<2:
                continue #no alternative
            outA_sets={}#transcripts supporting the different  junctions
            for tr,node_id in nA.suc.items(): 
                outA_sets.setdefault(node_id,set()).add(tr)
            unspliced=nA.end==self[junctions[0]].start
            alternative=([],outA_sets[junctions[0]]) if unspliced else (outA_sets[junctions[0]],[])
            logging.debug(f'checking node {i}: {nA} ({list(zip(junctions,[outA_sets[j] for j in junctions]))})')
            for idx,joi in enumerate(junctions[1:]): # start from second, as first does not have an alternative
                nB=self[joi]
                alternative=[{tr for tr in alternative[i] if self._pas[tr]>joi} for i in range(2)] #check that transcripts extend beyond nB
                logging.debug(alternative)
                weight=self.weights[:,list(outA_sets[joi])].sum(1)  
                found=[trL1.intersection(trL2) for trL1 in alternative for trL2 in inB_sets[joi]] #alternative transcript sets for the 4 types
                found.append(set.union(*alternative)-inB_sets[joi][0]-inB_sets[joi][1]) #5th type: mutually exclusive
                logging.debug(f'checking junction {joi} (tr={outA_sets[joi]}) and found {found} at B={inB_sets[joi]}')                
                for alt_type, alt in enumerate(found):
                    if alt:
                        alt_weight=self.weights[:,list(alt)].sum(1)
                        logging.debug(f'found loop {outA_sets[joi]} vs {alt} ({alt_type})')
                        yield weight, weight+alt_weight,nA.end,nB.start,alt_type
                alternative[0].update(outA_sets[joi]) #now transcripts supporting joi join the alternatives
                
                

                    
                


    def get_splice_coverage_old(self):
        'simpler version, but considering all transcripts spanning the junction (including transcripts where the splice sites are not exonic). Does not allow for classifying types'
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
