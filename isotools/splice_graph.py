from itertools import combinations
import pandas as pd
import numpy as np
from tqdm import tqdm
import logging
import scipy.stats as stats
import logging
from ._utils import pairwise, overlap
logger=logging.getLogger('isotools')

class SpliceGraph():
    def __init__(self, transcripts, strand):      
        self.strand=strand
        assert strand in '+-', 'strand must be either "+" or "-"'
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
                    if self.is_same_exon(trid,j1,j2) and self[j1].start<=exons[0][1] and self[j2].end >= exons[0][0]]
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


    def is_same_exon(self, tr_nr,j1,j2):
        #test if nodes j1 and j2 belong to same exon in transcript tr_nr
        for j in range(j1,j2):
            if tr_nr not in self[j].suc or self[j].suc[tr_nr]>j+1 or self[j].end!= self[j+1].start:
                return False
        return True
    
    def _count_introns(self,tr_nr,j1,j2):
        #count the number of junctions between j1 and j2
        logger.debug('counting introns of transcript %i between nodes %i and %i',tr_nr,j1,j2)
        delta=0
        if j1==j2:
            return 0
        assert tr_nr in self[j1].suc, f'transcript {tr_nr} does not contain node {j1}'
        while j1<j2:
            j_next=self[j1].suc[tr_nr]
            if j_next>j1+1 or self[j1].end != self[j1+1].start:
                delta+=1
            j1=j_next
        return delta

    def get_node_matrix(self):
        return np.array([[True if tss==j or trid in n.pre else False for j,n in enumerate(self)] for trid,tss in enumerate(self._tss)])    

    def find_fragments(self):
        truncated=set()
        contains={}
        starts=[set() for _ in range(len(self))]
        for trid,idx in enumerate(self._tss):
            starts[idx].add(trid)
        nodes=self.get_node_matrix()
        for trid,(tss,pas) in enumerate(zip(self._tss, self._pas)):
            if trid in truncated:
                continue
            contains[trid]={trid2 for trid2,(tss2,pas2) in enumerate(zip(self._tss, self._pas)) if trid2!=trid and tss2>= tss and pas2<= pas and all(nodes[trid2,tss2:pas2+1]==nodes[trid,tss2:pas2+1])}
            truncated.update(contains[trid])#those are not checked 
            #start_in=set() #transcripts that are contained until idx
            #contains[trid]=set()
            #while True:
            #    start_in.update(starts[idx])
            #    contains[trid].update({c for c in start_in if self._pas[c]==idx}) # add fully contained
            #    if idx == self._pas[trid]:
            #        contains[trid].remove(trid) #remove self
            #        truncated.update(contains[trid]) #those are not checked 
            #        break
            #    suc=self._graph[idx].suc
            #    start_in={c for c in start_in if c in suc and suc[c] ==suc[trid] } #remove transcripts that split at idx
            #    idx=suc[trid]
        fragments={}
        for big,smallL in contains.items():
            if big not in truncated:
                for trid in smallL:
                    delta1=self._count_introns(big,self._tss[big], self._tss[trid])
                    delta2=self._count_introns(big,self._pas[trid], self._pas[big])
                    fragments.setdefault(trid,[]).append((big,delta1,delta2) if self.strand=='+' else (big,delta2,delta1))
        return fragments


    def get_alternative_splicing(self, exons, alternative=None): 
        'compares exons to splice graph and returns list of novel splicing events'
        #returns a tuple
        # the sqanti category: 0=FSM,1=ISM,2=NIC,3=NNC,4=Novel gene
        # subcategories: a list of novel splicing events or splice_identical
        #t='novel exon','intron retention', 'novel exonic', 'novel intronic','splice identical','combination', 'fragment', 'extention'
        
        if alternative is not None and len(alternative)>0: # a list of tuples with (1) gene names and (2) junction numbers covered by other genes (e.g. readthrough fusion)
            category=4
            fusion_exons={int((i+1)/2) for j in alternative for i in j[1]}
            altsplice={'readthrough fusion':alternative}
        else:
            tr=self.search_transcript(exons)
            if tr:
                return 0,{'FSM':tr}        
            category=1
            altsplice={}
            fusion_exons=set()

        is_reverse=self.strand=='-'
        j1=next(j for j,n in enumerate(self) if n.end > exons[0][0])
        #j1: index of first segment ending after exon start (i.e. first overlapping segment)        
        try:
            j2=next(j-1 for j in range(j1,len(self)) if self[j].start >= exons[0][1])
        except StopIteration:
            j2=len(self)-1
        #j2: index of last segment starting befor exon end (i.e. last overlapping segment)   
     
        #check truncation at begining (e.g. low position)
        if (    len(exons)>1 and #no mono exon
                not any(j in self._tss for j in range(j1,j2+1)) and # no tss/pas within exon
                self[j1].start<= exons[0][0]): # start of first exon is exonic in ref
            j0=max(self._tss[trid] for trid in self[j1].pre)#j0 is the closest start node
            if any(self[j].end < self[j+1].start for j in range(j0,j1)): #assure there is an intron between closest tss/pas and exon
                end='3' if is_reverse else '5'
                altsplice.setdefault(f'{end}\' fragment',[]).append([self[j0].start, exons[0][0]]) #at start (lower position)


        for i,ex1 in enumerate(exons):                               
            ex2=None if i+1==len(exons) else exons[i+1]            
            if i not in fusion_exons:
                exon_altsplice,exon_cat=self._check_exon(j1,j2,i==0,is_reverse,ex1, ex2) #finds intron retention (NIC), novel exons, novel splice sites,  novel pas/tss  (NNC)
                category=max(exon_cat,category)
                for k,v in exon_altsplice.items(): 
                    altsplice.setdefault(k,[]).extend(v)
                # find j2: index of last segment starting befor exon2 end (i.e. last overlapping  segment)
                if ex2 is not None:
                    if j2+1<len(self):
                        j1,j2, junction_altsplice=self._check_junction(j2,ex1,ex2) #finds exon skipping and novel junction (NIC)
                        if junction_altsplice and i+1 not in fusion_exons:
                            category=max(2,category)
                            for k,v in junction_altsplice.items():
                                altsplice.setdefault(k,[]).extend(v)
                    else:
                        j1=len(self)

        
        #check truncation at end (e.g. high position)
        if ( len(exons)>1 and
                j2>=j1 and
                not any(j in self._pas for j in range(j1,j2+1)) and# no tss/pas within exon
                self[j2].end >= exons[-1][1]):# end of last exon is exonic in ref
            try:
                j3=min(self._pas[trid] for trid in self[j2].suc) #j3 is the next end node (pas/tss on fwd/rev)
            except ValueError:
                logger.error('\n'.join([str(exons), str(self._pas),str((j1,j2)),str([(j,n) for j,n in enumerate(self)])]))
                raise 
            if any(self[j].end < self[j+1].start for j in range(j2,j3)): #assure there is an intron between closest tss/pas and exon
                end='5' if is_reverse else '3'
                altsplice.setdefault(f'{end}\' fragment',[]).append([exons[-1][1],self[j3].end])
        
        if not altsplice:#all junctions are contained but not all in one transcript
            altsplice={'novel combination':[] }
            category=2


        return category,altsplice
            
    def _check_exon(self,j1,j2,is_first,is_reverse, e, e2=None):
        #checks whether exon is supported by splice graph between nodes j1 and j2
        #j1: index of first segment ending after exon start (i.e. first overlapping segment)
        #j2: index of last segment starting befor exon end (i.e. last overlapping  segment)

        logger.debug(f'check exon {e} between sg node {j1} and {j2}/{len(self)} (first={is_first},rev={is_reverse},e2={e2})')
        is_last=e2==None
        altsplice={}
        category=0
        if j1>j2: #e is not contained at all  -> novel exon (or TSS/PAS if first/last)
            category=3
            if is_first or is_last:
                altsplice={'novel PAS' if is_first==is_reverse else 'novel TSS':[e]}
            else:
                altsplice={'novel exon':[e]}
            j2=j1
        elif (is_first and is_last): # mono-exon            
            altsplice['mono-exon']=[]
            category=1
        else: # check splice sites
            if self[j1][0]!=e[0]: #first splice site missmatch
                if not is_first:
                    pos="intronic" if self[j1][0]>e[0] else "exonic"
                    kind='donor' if is_reverse else 'acceptor'
                    dist=e[0]-self[j1][0] #the distance to next junction
                    altsplice[f'novel {pos} splice {kind}']=[(e[0],dist)] 
                    category=3
                elif self[j1][0] > e[0] and not any(j in self._tss for j in range(j1,j2+1)): #exon start is intronic in ref
                    site='PAS' if is_reverse else 'TSS'
                    altsplice.setdefault(f'novel {site}',[]).append((e[0], self[j1][0]))
                    category=max(2,category)  
            if self[j2][1]!=e[1]:#second splice site missmatch
                if not is_last:
                    pos="intronic" if self[j2][1]<e[1] else "exonic"
                    kind='acceptor' if is_reverse else 'donor'
                    dist= e[1]-self[j2][1]
                    altsplice.setdefault(f'novel {pos} splice {kind}',[]).append((e[1],dist))
                    category=3
                elif self[j2][1] < e[1] and not any(j in self._pas for j in range(j1,j2+1)): #exon end is intronic in ref
                    site='TSS' if is_reverse else 'PAS'
                    altsplice.setdefault(f'novel {site}',[]).append((self[j2][1],e[1]))
                    category=max(2,category)    

        #find intron retentions
        if j1<j2:
            ret_introns=[]
            troi=set(self[j1].suc.keys()).intersection(self[j2].pre.keys())
            if troi:
                j=j1
                while j<j2:
                    nextj=min(js for trid,js in self[j].suc.items() if trid in troi)
                    if self[nextj].start-self[j].end>0:
                        ret_introns.append((self[nextj].start,self[j].end))
                    j=nextj
                if ret_introns:                
                    altsplice['intron retention']=ret_introns
                    category=max(2,category)
        logger.debug(f'check exon {e} resulted in {altsplice}')
        return altsplice, category

    def _check_junction(self, j2, e, e2):
        #1) check presence e1-e2 junction in ref (-> if not exon skipping or novel junction)
        #2) find j3: first node overlapping e2 and j4: last node overlapping e2
        altsplice={}
        introns=[]
        for j3 in range(j2+1,len(self)-1): 
            if j3 in self[j3-1].suc.values() and self[j3-1][1]<self[j3][0]: #"real" junction 
                introns.append([self[j3-1][1],self[j3][0]])
            if self[j3][1] > e2[0]:
                break
        else:
            j3=len(self) #e2 does not overlap splice graph
        if j3<len(self) and j3 not in self[j2].suc.values(): #direkt e1-e2 junction not present
            if set.intersection(set(self[j2].suc.keys()),set(self[j3].pre.keys())) and len(introns)>1: #path from e1 to e2 present
                #assert len(introns)>1,'logical issue - exon skipping should involve two introns'
                for i1,i2 in pairwise(introns):
                    altsplice.setdefault('exon skipping',[]).append([i1[1], i2[0]])
            elif e[1]==self[j2].end and e2[0]==self[j3].start: #e1-e2 path is not present, but splice sites are
                altsplice.setdefault('novel junction',[]).append([e[1],e2[0]]) #for example mutually exclusive exons spliced togeter
        
        logger.debug(f'check junction {e[1]} - {e2[0]} resulted in {altsplice}')
        #find last node overlapping e2
        try:
            j4=next(j-1 for j in range(j3,len(self)) if self[j].start >= e2[1])
        except StopIteration:
            j4=len(self)-1
        
        return j3,j4, altsplice
            
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

        
    def find_splice_sites(self, exons):
        j=0
        sites=np.zeros((len(exons)-1)*2, dtype=bool)
        for i,(e1,e2) in enumerate(pairwise(exons)):
            while(self[j].end<e1[1]):
                j+=1
                if j==len(self):
                   return sites
            if self[j].end==e1[1]:
                sites[i*2]=True
            while(self[j].start<e2[0]):
                j+=1
                if j==len(self):
                   return sites
            if self[j].start==e2[0]:
                sites[i*2+1]=True
        return sites

    def get_overlap(self,exons):
        ol=0
        j=0
        for e in exons:
            while(self[j].end<e[0]):#no overlap
                j+=1
                if j==len(self):
                   return ol
            while(self[j].start<e[1]):
                i_end = min(e[1], self[j].end)
                i_start = max(e[0], self[j].start)
                ol += (i_end-i_start)
                j+=1
                if j==len(self):
                   return ol
        return ol

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


    def find_ts_candidates(self, coverage):        
        for i, gnode in enumerate(self._graph[:-1]):
            if self._graph[i+1].start==gnode.end: #jump candidates: introns that start within an exon
                jumps={idx:n for idx,n in gnode.suc.items() if n>i+1 and self._graph[n].start==self._graph[n-1].end}
                #find jumps (n>i+1) and check wether they end within an exon begin(jumptarget)==end(node before)
                jump_weight={}
                for idx, target in jumps.items():
                    jump_weight.setdefault(target, [0,[]])
                    jump_weight[target][0]+=coverage[:,idx].sum(0)
                    jump_weight[target][1].append(idx)

                for target, (w,idx) in jump_weight.items():
                    long_idx=set(idx for idx,n in gnode.suc.items() if n==i+1) & set(idx for idx,n in self[target].pre.items() if n==target-1)
                    try:
                        longer_weight=coverage[:,list(long_idx)].sum()
                    except IndexError:
                        print(long_idx)
                        raise
                    yield gnode.end, self[target].start,w, longer_weight,  idx
        
    def splice_dependence(self, sidx,min_cov,coverage, ignore_unspliced=True):
        'experimental'
        # TODO - outdated - remove when replaced by splice bubbles'
        tr_cov=coverage[sidx,:].sum(0)
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
            if ignore_unspliced:
                nojump_ids= {trid for trid in jstart.intersection(jend) if self._is_spliced(trid,j[0], j[1])}
            else:
                nojump_ids=jstart.intersection(jend)
            if tr_cov[list(nojump_ids)].sum()>=min_cov:
                nonjumps[j]=nojump_ids
        for j1 in nonjumps:
            for j2 in (j for j in nonjumps if j[0]>j1[1]):
                double_negative=list(nonjumps[j1].intersection(nonjumps[j2]))
                #if ignore_unspliced:
                #    double_negative=[trid for trid in double_negative if self._is_spliced(trid,j1[0], j2[1])]
                ma=np.array([
                    tr_cov[list(jumps[j1].intersection(jumps[j2]))].sum(),
                    tr_cov[list(jumps[j1].intersection(nonjumps[j2]))].sum(),
                    tr_cov[list(nonjumps[j1].intersection(jumps[j2]))].sum(),
                    tr_cov[double_negative].sum()]).reshape((2,2))
                if any(ma.sum(0)<min_cov) or any(ma.sum(1)<min_cov):
                    continue
                #oddsratio, pvalue = fisher_exact(ma)                
                #log.debug('j1={}, j2={}, ma={}, or={}, p={}'.format(j1,j2,list(ma),oddsratio, pvalue))
                yield ma,(self[j1[0]].end, self[j1[1]].start), (self[j2[0]].end, self[j2[1]].start)

    def _is_spliced(self,trid,ni1,ni2):
        'checks if transcript is spliced (e.g. has an intron) between nodes ni1 and ni2'
        if any(self[i].end<self[i+1].start for i in range(ni1, ni2)): #all transcripts are spliced
            return True
        if all(trid in self[i].suc for i in range(ni1, ni2)): 
            return False
        return True

    def _get_next_spliced(self, trid, node):
        'find the next spliced node for given transcript'
        while node!=self._pas[trid]:
            next_node=self[node].suc[trid] #raises error if trid not in node
            if self[next_node].start > self[node].end:
                return next_node
            node=next_node
        return None

    def _get_exon_end(self, trid, node):
        'find the end of the exon to which node belongs for given transcript'
        while node!=self._pas[trid]:
            next_node=self[node].suc[trid] #raises error if trid not in node
            if self[next_node].start > self[node].end:
                return node
            node=next_node
        return node
    def _get_exon_start(self, trid, node):
        'find the start of the exon to which node belongs for given transcript'
        while node!=self._tss[trid]:
            next_node=self[node].pre[trid] #raises error if trid not in node
            if self[next_node].end < self[node].start:
                return node
            node=next_node
        return node
    

    def find_splice_bubbles(self, include_starts=True):
        '''search for alternative paths in the splice graphs ("bubbles"), 
        e.g. combinations of nodes A and B with more than one path from A to B.
        returns the transcript numbers and node ids of most direkt paths and alternative paths respectivly and type of alternative'''
         # alternative types: intron retention, alternative splice site at left and right, exon skipping, mutually exclusive
        alt_types=['ES','5AS', '3AS','IR','ME'] if self.strand=='-' else ['ES','3AS', '5AS','IR','ME']
        inB_sets=[(set(),set())] #list of spliced and unspliced transcripts joining in B
        #node_matrix=self.get_node_matrix()
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
            #nC_dict aims to avoid recalculation of nC for ME events
            nC_dict={} # tr -> node at start of 2nd exon C for tr such that there is one exon (B) (and both flanking introns) between nA and C; None if transcript ends
            me_alt_seen=set()
            logger.debug(f'checking node {i}: {nA} ({list(zip(junctions,[outA_sets[j] for j in junctions]))})')
            for j_idx,joi in enumerate(junctions[1:]): # start from second, as first does not have an alternative
                alternative=[{tr for tr in alternative[i] if self._pas[tr]>joi} for i in range(2)] #check that transcripts extend beyond nB
                logger.debug(alternative)
                found=[trL1.intersection(trL2) for trL1 in alternative for trL2 in inB_sets[joi]] #alternative transcript sets for the 4 types
                #found.append(set.union(*alternative)-inB_sets[joi][0]-inB_sets[joi][1]) #5th type: mutually exclusive (outdated handling of ME for reference)
                logger.debug(f'checking junction {joi} (tr={outA_sets[joi]}) and found {found} at B={inB_sets[joi]}')                
                for alt_tid, alt in enumerate(found):
                    if alt:
                        yield list(outA_sets[joi]),list(alt),i,joi,alt_types[alt_tid]
                #me_alt=set.union(*alternative)-inB_sets[joi][0]-inB_sets[joi][1] #search 5th type: mutually exclusive   
                me_alt=alternative[0]-inB_sets[joi][0]-inB_sets[joi][1] #search 5th type: mutually exclusive - needs to be spliced
                if me_alt-me_alt_seen: #there is at least one novel alternative transcript
                    #for ME we need to find (potentially more than one) node C where the alternatives rejoin
                    for tr in me_alt:
                        if tr not in nC_dict:
                            nC_dict[tr]=self._get_next_spliced(tr,nA.suc[tr])
                            if nC_dict[tr] is None: #transcript end, no node C
                                me_alt_seen.add(tr) #those are not of interest for ME
                    inC_sets={}
                    for nB_i in junctions[j_idx+1:]:
                        for tr in outA_sets[nB_i]:
                            if tr not in nC_dict:
                                nC_dict[tr]=self._get_next_spliced(tr,nB_i)
                            nC_i=nC_dict[tr]
                            if nC_i is not None:
                                if nB_i==joi:
                                    inC_sets.setdefault(nC_i,set()).add(tr)
                                elif nC_i in inC_sets:
                                    inC_sets[nC_i].add(tr)
                        if not inC_sets:
                            break
                    for nC_i, me_prim in sorted(inC_sets.items()):
                        found_alt={tr for tr in me_alt if nC_dict[tr]==nC_i}
                        if found_alt-me_alt_seen:
                            yield (list(me_prim), list(found_alt),i,nC_i,'ME')
                            me_alt=me_alt-found_alt
                            me_alt_seen.update(found_alt)                    
                alternative[0].update(outA_sets[joi]) #now transcripts supporting joi join the alternatives
        if include_starts:
            yield from self.find_alternative_starts()
        
    def find_alternative_starts(self):
        #analog to splice bubbles
        #actually finds first/last splice site
        tss={}
        pas={}
        tss_start={}
        pas_end={}
        for tr,(start1, end1) in enumerate(zip(self._tss, self._pas)):
            start2=self._get_exon_end(tr,start1)
            if start2 != end1:
                tss.setdefault(start2,set()).add(tr)
                tss_start[start2]=min(start1, tss_start.get(start2,start1))
            end2=self._get_exon_start(tr,end1)
            if end2 != start1:
                pas.setdefault(end2,set()).add(tr)
                pas_end[end2]=max(end1, pas_end.get(end2,end1))
        for n,tr_set in tss.items(): # TODO: maybe also add compatible alternatives: end after tss /start before pas
            alt_tr=[tr for tr in range(len(self._tss)) if tr not in tr_set and self._pas[tr] > n]
            yield (alt_tr,list(tr_set),tss_start[n],n,'TSS' if self.strand=='+' else 'PAS')
        for n,tr_set in pas.items():
            alt_tr=[tr for tr in range(len(self._tss)) if tr not in tr_set and self._tss[tr] < n]
            yield (alt_tr,list(tr_set),pas_end[n],n,'PAS' if self.strand=='+' else 'TSS')


    def _get_all_exons(self, nodeX, nodeY, tr):    
        'get all exons from nodeX to nodeY for transcripts tr'
        #TODO: add option to extend first and last exons
        n=max(nodeX, self._tss[tr]) # if tss>nodeX start there
        if tr not in self[n].pre and self._tss[tr]!=n: #nodeX is not part of tr
            for i in range(n,nodeY+1):
                if tr in self[n].suc:
                    n=i
                    break
            else: 
                return []
        if n>nodeY:
            return []
        exons=[[self[n].start, self[n].end]]
        while n<nodeY:
            try:
                n=self[n].suc[tr]
            except KeyError: #end of transcript before nodeY was reached
                break
            if self[n].start==exons[-1][1]:
                exons[-1][1]=self[n].end
            else:
                exons.append([self[n].start, self[n].end])
        return [tuple(e) for e in exons]
        
    def __getitem__(self, key):
        return self._graph[key]
    
    def __len__(self):
        return len(self._graph)

      

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
