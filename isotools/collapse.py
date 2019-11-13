

import isotools.utils
import isotools.splice_graph

def collapse_transcript_of_gene(gene, fuzzy_junction, repair_5truncation, repair_3trunkation, rename=False):    
    collapsed=dict()
    if fuzzy_junction>0:
        try:
            isotools.splice_graph.unify_junctions_of_gene(gene, fuzzy_junction)
        except Exception as e:
            print(e)
            print(gene)
            raise
    for i,(tr1_name, tr1) in enumerate(gene.data['transcripts'].items()):
        for tr2_name, tr2 in collapsed.items():
            exons1=tr1['exons']
            exons2=tr2['exons']
            trunkated=is_truncation(exons1, exons2)
            if trunkated:
                #trunkation detected, find out the exonic distance
                relation, first, last=trunkated
                if first==0:
                    delta1=sum( [exons2[i][1] -exons2[i][0] for i in range(relation[0][0][0])])
                else:
                    delta1=-sum( [exons1[i][1] -exons1[i][0] for i in range(first)])
                delta1+=exons1[first][0]-exons2[relation[first][0][0]][0]
                if last==len(exons1)-1:
                    delta2=sum( [exons2[i][1] -exons2[i][0] for i in range(relation[last][0][0]+1, len(exons2))])
                else:
                    delta2=-sum( [exons1[i][1] -exons1[i][0] for i in range(last+1, len(exons1))])
                delta2-=exons1[last][1]-exons2[relation[last][0][0]][1]
                #todo: check thresholds
                collapsed[tr2_name]=merge_transcripts(tr1,tr2)
                break                        
        else:
            collapsed[tr1_name]=tr1
    gene.data['transcripts']=collapsed


def merge_transcripts(tr1, tr2, new_name=None):
    #print('merging {} {}'.format(tr1, tr2))
    merged={'exons':list()}
    e2iter=iter(tr2['exons'])
    e1enum=enumerate(tr1['exons'])
    e2=next(e2iter)
    for i,e1 in e1enum:
        if isotools.utils.overlap(e1,e2):
            merged['exons'].append( ( min(e1[0], e2[0]), max(e1[1], e2[1]) ) )
        elif e1[0]<e2[0]:
            merged['exons'].append(e1)
        else:
            merged['exons'].append(e2)
            try:
                e2=next(e2iter)
            except StopIteration:
                merged['exons']+=[e for j,e in e1enum if j> i]
                break
    else:
        merged['exons'].append(e2)
        merged['exons']+=list(e2iter)

    for n in tr1.keys()|tr2.keys():
        if n not in ['exons', 'support']:
            merged[n]=tr1.get(n,[])+tr2.get(n,[]) #carefull, might not be possible for all datatypes

    return merged
        




def is_truncation(exons1, exons2, debug=False):    
    relation=isotools.utils.get_relation(exons1, exons2)
    if any(len(r) > 1 for r in relation):#
        if debug: print('one exon from 1 has several corresponding exons from 2')
        return False  
    if not any(r for r in relation ): #
        if debug: print('no relation')
        return False 
    first=next(iter([i for i,r in enumerate(relation) if r]))
    last=len(relation)-next(iter([i for i,r in enumerate(reversed(relation)) if r]))-1
    if sum(1 for r in relation if r) != last-first+1:  # 
        if debug: print('some exons in 1 do not have corresponding exons in 2')
        return False
    if first>0 and relation[first][0][0]>0:
        if debug: print('start exons do not correspond to each other')
        return False # 
    if last < len(exons1)-1 and relation[last][0][0]<len(exons2)-1 :
        if debug: print('last exons do not correspond to each other')
        return False # 
    if last!=first and ( #more than one exon
            (relation[first][0][1] & 1 ==0) or # 2nd splice junction must match
            (relation[last][0][1] & 2 == 0)): # fst splice junction must match
        if debug: print('more than one exon, 2nd splice junction must match, fst splice junction must match')
        return False        
    if (relation[first][0][0]!=first and # one of the transcripts has extra exons at start
            (first==0 ) == (exons1[first][0] < exons2[relation[first][0][0]][0])): # check the begin of the shorter
        if debug: print('first exon of trunkated version is actually longer')
        return False
    if (relation[last][0][0] != len(exons2)-1 and # one of the transcripts has extra exons at the end
            (last==len(exons1)-1) == (exons1[last][1] > exons2[relation[last][0][0]][1])):  #check the end of the shorter
        if debug: print('last exon of trunkated version is actually longer')
        return False
    if last-first > 3 and any(relation[i][0][1]!=3 for i in range(first+1, last)): #
        if debug: print('intermediate exons do not fit')
        return False
    if debug: print('all filters passed')
    return relation, first, last



def merge_interval(e1, e2, must_overlap=False):
    if must_overlap and not isotools.utils.overlap(e1, e2):
        raise ValueError('intervals {} and {} do not overlap'.format(e1, e2))
    return (min(e1[0], e2[0]), max(e1[1], e2[1]))

def shared_exons(tr1, tr2):
    rel = isotools.utils.get_relation(tr1, tr2)
    same = {r[0] for x in rel for r in x if r[1] == 3}

    return len(same)

def splice_identical(tr1, tr2):
    if len(tr1) != len(tr2):  # different number of exons
        return False
    if len(tr1) == 1 and isotools.utils.overlap(tr1[0], tr2[0]):  # single exon genes
        return True
    if tr1[0][1] != tr2[0][1] or tr1[-1][0] != tr2[-1][0]:  # check first and last exons
        return False
    for e1, e2 in zip(tr1[1:-1], tr2[1:-1]):  # check other exons
        if e1[0] != e2[0] or e1[1] != e2[1]:
            return False
    return True
