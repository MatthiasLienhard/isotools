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


def import_gtf_transcripts(refseq_fn, genes=None):
    # returns a dict interval trees for the genes, each containing the splice graph
    gff = TabixFile(refseq_fn)
    chrom_ids = get_gff_chrom_dict(gff)
    exons = dict()  # transcript id -> exons
    transcripts = dict()  # gene_id -> transcripts
    skipped = set()
    if genes is None:
        genes=dict()
    #takes quite some time... add a progress bar?
    for line in tqdm(gff.fetch(), smoothing=.1):  # parser=pysam.asGTF()): parser is strange and not helpfull
        ls = line.split(sep="\t")
        if ls[0] not in chrom_ids:
            #warnings.warn('unknown chromosome id :'+ls[0])
            # this happens at all NW/NT scaffolds
            continue
        chrom = chrom_ids[ls[0]]
        genes.setdefault(chrom, IntervalTree())
        info = dict([pair.split("=", 1) for pair in ls[8].split(";")])
        start, end = [int(i) for i in ls[3:5]]
        start -= 1  # to make 0 based
        if 'Name' not in info:
            for alt in 'standard_name', 'gene', 'product':
                try:
                    info['Name'] = info[alt]
                    break
                except KeyError:
                    pass
        if ls[2] == "exon":
            # print(line)
            if "Parent" in info.keys():
                gtf_id = info['Parent']
                exons.setdefault(gtf_id, list()).append((start, end))
            else:  # should not happen if GTF is OK
                warnings.warn("Skipping exon without id: "+line)
        # elif ls[2] in ['gene', 'pseudogene']:
        elif ls[2] == 'gene':
            info['strand'] = ls[6]
            info['chromosome'] = chrom
            genes[chrom][start:end] = info
        # those denote transcripts
        elif all([v in info for v in ['Parent', "ID", 'Name']]) and info['Parent'].startswith('gene'):
            transcripts.setdefault(info["Parent"], list()).append(
                (info["Name"], info["ID"]))
        else:
            skipped.add(ls[2])
    #print('skipped the following categories: {}'.format(skipped))
    # sort the exons
    for tid in exons.keys():
        exons[tid].sort()
    # add transcripts to genes
    for chrom in genes:
        for gene in genes[chrom]:
            g_id = gene.data['ID']
            t_ids = transcripts.get(g_id, [(gene.data['Name'], g_id)])
            for t_name, t_id in t_ids:
                try:
                    gene.data.setdefault('transcripts', dict())[t_name] = {'exons':exons[t_id]}
                except KeyError:
                    # genes without transcripts get a single exons transcript
                    gene.data['transcripts'] = {t_name: {'exons':[tuple(gene[:2])]}}
                    
    return genes

def import_bam_transcripts(bam_fn,gene_prefix='PB_',  genes=None, n_lines=None):
    genes = dict()
    gene_number=0
    merged_number=0
    with AlignmentFile(bam_fn, "rb") as align:
        
        stats = align.get_index_statistics()
        # try catch if sam/ no index
        total_reads = sum([s.mapped for s in stats])
        if n_lines is not None:
            total_reads=n_lines
        for i, read in tqdm(enumerate(align.fetch()), total=total_reads, unit='transcripts'):
            chrom = read.reference_name
            strand = '-' if read.is_reverse else '+'
            genes.setdefault(chrom, IntervalTree())
            exons = junctions_from_read(read)
            pos = exons[0][0], exons[-1][1]
            found=None
            for gene in genes[chrom].overlap(*pos):                
                if gene.data['strand'] == strand and any([is_same_gene(exons, tr['exons']) 
                            for tr in gene.data['transcripts'].values()]):
                    if found is None:
                        tr_name='{}_{}'.format(gene.data['Name'],len(gene.data['transcripts'])+1)
                        gene.data['transcripts'][tr_name]={'exons':exons, 'source':[read.query_name]}
                        if pos[0]<gene.begin or pos[1]>gene.end:
                            merged_gene=Interval(min(pos[0],gene.begin), max(pos[1], gene.end), gene.data)
                            genes[chrom].add(merged_gene)
                            genes[chrom].remove(gene)    
                            found=merged_gene
                        else:
                            found=gene
                    else:
                        #merge (current) gene and (already) found_gene, as both overlap with new transcript
                        #print('this happens!!')
                        merged_number+=1
                        found.data['transcripts'].update(gene.data['transcripts'])
                        if gene.begin<found.begin or gene.end > found.end:
                            merged_gene=Interval(min(found.begin,gene.begin), max(found.end, gene.end), found.data)
                            genes[chrom].add(merged_gene)
                            genes[chrom].remove(found)    
                            found=merged_gene
                        genes[chrom].remove(gene)                            
            if found is None: #new gene
                gene_number+=1
                gname=gene_prefix+str(gene_number)
                info={'strand':strand, 'Name':gname,'source':{gname+'_1':[read.query_name]}, 
                        'transcripts':{gname+'_1':{'exons':exons,'source':[read.query_name]}}}
                genes[chrom].addi(*pos, info)
            if n_lines==i:
                break
        print('imported {} transcripts in {} genes'.format(total_reads,gene_number-merged_number ))
    return genes

def add_splice_graphs(genes, force=False):
    with tqdm(total=sum([len(tree) for tree in genes.values()]), unit='genes') as pbar:     
        for chrom, tree in genes.items():
            pbar.set_postfix(chr=chrom)
            for gene in tree:
                pbar.update(1)
                if not force and 'splice_graph' in gene.data:
                    continue
                gene.data['splice_graph']=get_splice_graph(gene)
                                            


                
def get_splice_graph(gene):
    return(SpliceGraph((tr['exons'] for tr in gene.data['transcripts'].values()), gene.end))

def get_read_sequence(bam_fn,reg=None, name=None):
    with AlignmentFile(bam_fn, "rb") as align:
        for i, read in enumerate(align.fetch(region=reg)):
            chrom = read.reference_name
            strand = '-' if read.is_reverse else '+'
            exons = junctions_from_read(read)
            if name is None or name == read.query_name:
                yield  read.query_name, read.query_sequence, chrom,strand, exons


def collapse_transcripts(genes,  fuzzy_junction=5, repair_5truncation=10000, repair_3trunkation=100, rename=True):
    with tqdm(total=sum([len(tree)for tree in genes.values()]), unit='genes') as pbar:     
        for chrom, tree in genes.items():
            for gene in tree:
                collapse_transcript_of_gene(gene, fuzzy_junction, repair_5truncation, repair_3trunkation, rename)
                pbar.update(1)

def merge_transcripts(tr1, tr2, new_name=None):
    #print('merging {} {}'.format(tr1, tr2))
    merged={'exons':list()}
    e2iter=iter(tr2['exons'])
    e1enum=enumerate(tr1['exons'])
    e2=next(e2iter)
    for i,e1 in e1enum:
        if overlap(e1,e2):
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
    relation=get_relation(exons1, exons2)
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



def collapse_transcript_of_gene(gene, fuzzy_junction, repair_5truncation, repair_3trunkation, rename=False):    
    collapsed=dict()
    if fuzzy_junction>0:
        try:
            unify_junctions_of_gene(gene, fuzzy_junction)
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

def unify_junctions_of_gene(gene, fuzzy_junction=5):
    if 'splice_graph' in gene.data:
        splice_graph=gene.data['splice_graph']
    else:
        splice_graph=get_splice_graph(gene)
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
    for tr in gene.data['transcripts'].values():
        tr['exons']=[(new_start[e[0]], new_end[e[1]]) for e in tr['exons']]
    

            



def add_reference_support(genes, ref_genes):
    n = sum([len(gene.data['transcripts'])
                for tree in genes.values() for gene in tree])
    with tqdm(total=n, unit='transcripts') as pbar:
        for chrom, tree in genes.items():
            pbar.set_postfix(chr=chrom)
            for gene in tree:
                for tr in gene.data['transcripts'].values():
                    tr['support']=get_support(tr['exons'], ref_genes, chrom=chrom,is_reverse=gene.data['strand'] == '-')
                    pbar.update(1)


#utils
def num_transcripts(genes):
    return sum([len(gene.data['transcripts']) for tree in genes.values() for gene in tree])

def get_support(exons, ref_genes, chrom, is_reverse):
    support = {'ref_gene': 'NA', 'ref_transcript': 'NA', 'ref_tss': 0, 'ref_pas': 0,
                        'ref_len': 0, 'ref_nSJ': 0, 'exI': 0, 'sjI': 0, 'exIoU': 0, 'sjIoU': 0, 'sType': ['novel/unknown']}
    if chrom not in ref_genes:
        return support
    ref_genes_ol = ref_genes[chrom][exons[0][0]: exons[-1][1]]
    support_dict = compute_support(ref_genes_ol, exons)
    try:
        # best_idx=np.argmax(zip(support_dict['sjIoU'],support_dict['exIoU']))
        best_idx = max(enumerate(
            zip(support_dict['sjIoU'], support_dict['exIoU'])), key=lambda x: x[1])[0]
        # https://stackoverflow.com/questions/2474015/getting-the-index-of-the-returned-max-or-min-item-using-max-min-on-a-list
    except ValueError:
        return support
    else:
        for n, v in support_dict.items():
            support[n]=v[best_idx]
        support["sType"]=get_alternative_splicing(support_dict, best_idx,exons, ref_genes_ol, is_reverse)
    return support

def get_alternative_splicing(support_dict, best_idx,exons, ref_genes_ol, is_reverse):
    if support_dict["sjIoU"][best_idx] == 1:
        splice_type = ['splice_identical']
    elif support_dict["exIoU"][best_idx] == 0:
        splice_type = ['novel/unknown']
    else:
        best_tid = support_dict['ref_transcript'][best_idx]
        for g in ref_genes_ol:
            if best_tid in g.data['transcripts']:
                ref_exons = g.data['transcripts'][best_tid]['exons']
                splice_type = get_splice_type(
                    ref_exons, exons, is_reverse)
                break
        else:
            warnings.warn(
                'cannot find the transcript -- this should never happen')
        covered_gene_ids = {id for id, iou in zip(
            support_dict['ref_gene'], support_dict['sjIoU']) if iou > 0}
        covered_genes = {
            g for g in ref_genes_ol if g.data['Name'] in covered_gene_ids}
        if len(covered_genes) > 1:
            # remove genes contained in other genes
            contained = []
            for g1, g2 in combinations(covered_genes, 2):
                if g1[0] <= g2[0] and g1[1] >= g2[1]:
                    # g2 is contained in g1
                    contained.append(g2)
                elif g1[0] >= g2[0] and g1[1] <= g2[1]:
                    # g1 is contained in g2
                    contained.append(g1)
            for g1, g2 in combinations([g for g in covered_genes if g not in contained], 2):
                if not overlap(g1[:2], g2[:2]):
                    splice_type.setdefault('fusion_gene', []).append(
                        '{}|{}'.format(g1.data['Name'], g2.data['Name']))
    return splice_type
    #support["sType"].append(splice_type)
    #else:
    #    for n, v in support_default.items():
    #        support.setdefault(n, []).append(v)
    #return support

def is_same_gene(tr1, tr2, spj_iou_th=0, reg_iou_th=.9):
    # current definition of "same gene": at least one shared splice site
    # or at least 90% exonic overlap
    spj_i, reg_i=get_intersects(tr1, tr2)
    total_spj=(len(tr1)+len(tr2)-2)*2
    total_len=sum([e[1]-e[0] for e in tr2+tr1])
    spj_iou=spj_i/(total_spj-spj_i) if total_spj>0 else 0
    reg_iou=reg_i/(total_len-reg_i)
    if spj_iou>spj_iou_th or reg_iou>reg_iou_th:
        return True
    return False

def get_gff_chrom_dict(gff):
    # fetch chromosome ids
    chrom = {}
    for c in gff.contigs:
        #print ("---"+c)
        for line in gff.fetch(c, 1, 2):
            if line[1] == "C":
                ls = line.split(sep="\t")
                if ls[2] == "region":
                    info = dict([pair.split("=")
                                    for pair in ls[8].split(";")])
                    if "chromosome" in info.keys():
                        chrom[ls[0]] = info["chromosome"]
                        # print(line)
    return(chrom)

'''
def add_exon_to_splice_graph(splice_graph, exon, predecessor, successor):
    ol_exons = splice_graph.overlap(*exon)
    if not ol_exons:  # no overlapping exons so far, just add to splice tree
        splice_graph.add(Interval(*exon, (predecessor, successor)))
        return
    links = {exon[0]: predecessor, exon[1]: successor}
    # the keys contain the positions, where exons start or end
    # the values contain the end/start positions of predecessors and successors (as a set)
    # for exon splits, key equals one of the values
    for ol in ol_exons:
        if ol.begin in links:
            links[ol.begin].update(ol.data[0])
        else:
            links[ol.begin] = ol.data[0]
            if exon[0] < ol.begin and ol.begin < exon[1]:
                links[ol.begin].add(ol.begin)
            else:
                links[exon[0]].add(exon[0])
        if ol.end in links:
            links[ol.end].update(ol.data[1])
        else:
            links[ol.end] = ol.data[1]
            if exon[0] < ol.end and ol.end <= exon[1]:
                links[ol.end].add(ol.end)
            else:
                links[exon[1]].add(ol.end)
        splice_graph.remove(ol)
    links_list = sorted(links.items())
    for i, end in enumerate(links_list[1:]):
        start = links_list[i]
        splice_graph.add(Interval(start[0], end[0],
                                    ({s_pos for s_pos in start[1] if s_pos <= start[0]},
                                    {e_pos for e_pos in end[1] if e_pos >= end[0]})))


def check_isoseq_bam(ref_genes,bam_fn='/project/42/pacbio/hecatos_isoseq/05-gmap/all_merged_collapsed_isoforms.bam'):
    support = {k: list() for k in ['name', 'chrom', 'tss', 'pas', 'len', 'nSJ', 'ref_gene',
                                    'ref_transcript', 'ref_tss', 'ref_pas', 'ref_len', 'ref_nSJ', 'exI', 'sjI', 'sType']}
    with AlignmentFile(bam_fn, "rb") as align:
        stats = align.get_index_statistics()
        # try catch if sam/ no index
        total_reads = sum([s.mapped for s in stats])
        for read in tqdm(align.fetch(), total=total_reads, unit='transcripts'):
            # continue if not mapped??
            add_support(ref_genes, support,
                                chrom=read.reference_name,
                                name=read.query_name[:read.query_name.find(
                                    '|')],
                                exons=junctions_from_read(read),
                                is_reverse=read.is_reverse)
    return support
'''
def add_support(ref_genes, support, chrom, name, exons, is_reverse=False):
    support_default = {'ref_gene': 'NA', 'ref_transcript': 'NA', 'ref_tss': 0, 'ref_pas': 0,
                        'ref_len': 0, 'ref_nSJ': 0, 'exI': 0, 'sjI': 0, 'exIoU': 0, 'sjIoU': 0, 'sType': ['novel/unknown']}
    support['name'].append(name)
    support['chrom'].append(chrom)
    support['tss'].append(
        exons[0][0] if not is_reverse else exons[-1][1])
    support['pas'].append(
        exons[0][0] if is_reverse else exons[-1][1])
    support['len'].append(sum([e[1]-e[0] for e in exons]))
    support['nSJ'].append(len(exons)*2-2)

    if chrom in ref_genes:
        ref_genes_ol = ref_genes[chrom][exons[0][0]: exons[-1][1]]
        support_dict = compute_support(ref_genes_ol, exons)
        try:
            # best_idx=np.argmax(zip(support_dict['sjIoU'],support_dict['exIoU']))
            best_idx = max(enumerate(
                zip(support_dict['sjIoU'], support_dict['exIoU'])), key=lambda x: x[1])[0]
            # https://stackoverflow.com/questions/2474015/getting-the-index-of-the-returned-max-or-min-item-using-max-min-on-a-list
        except ValueError:
            for n, v in support_default.items():
                support.setdefault(n, []).append(v)
        else:
            for n, v in support_dict.items():
                support.setdefault(n, []).append(v[best_idx])
            if support_dict["sjIoU"][best_idx] == 1:
                splice_type = ['splice_identical']
            elif support_dict["exIoU"][best_idx] == 0:
                splice_type = ['novel/unknown']
            else:
                best_tid = support_dict['ref_transcript'][best_idx]
                for g in ref_genes_ol:
                    if best_tid in g.data['transcripts']:
                        ref_exons = g.data['transcripts'][best_tid]['exons']
                        splice_type = get_splice_type(
                            ref_exons, exons, is_reverse)
                        break
                else:
                    warnings.warn(
                        'cannot find the transcript -- this should never happen')
                covered_gene_ids = {id for id, iou in zip(
                    support_dict['ref_gene'], support_dict['sjIoU']) if iou > 0}
                covered_genes = {
                    g for g in ref_genes_ol if g.data['Name'] in covered_gene_ids}
                if len(covered_genes) > 1:
                    # remove genes contained in other genes
                    contained = []
                    for g1, g2 in combinations(covered_genes, 2):
                        if g1[0] <= g2[0] and g1[1] >= g2[1]:
                            # g2 is contained in g1
                            contained.append(g2)
                        elif g1[0] >= g2[0] and g1[1] <= g2[1]:
                            # g1 is contained in g2
                            contained.append(g1)
                    for g1, g2 in combinations([g for g in covered_genes if g not in contained], 2):
                        if not overlap(g1[:2], g2[:2]):
                            splice_type.setdefault('fusion_gene', []).append(
                                '{}|{}'.format(g1.data['Name'], g2.data['Name']))
            support["sType"].append(splice_type)
    else:
        for n, v in support_default.items():
            support.setdefault(n, []).append(v)

def alt_splice_fraction(support=None):
    stype = dict(zip(support['name'], support['sType']))
    type_list = list(itertools.chain.from_iterable(stype.values()))
    type_counts = sorted(
        list(zip(*np.unique(type_list, return_counts=True))), key=lambda x: x[1])
    total = len(stype)
    type_fraction = {n: round(v/total*100, 2) for n, v in type_counts}
    return(type_fraction)

def barplot(data, ax=None):
    if ax is None:
        ax = plt.gca()
    ax.barh(range(len(data)), list(
        data.values()), align='center')

    ax.set_yticks(range(len(data)))
    ax.set_yticklabels(list(data.keys()))
    return ax

    # out_fn=out_fn=os.path.splitext(fn)[0]+'_altsplice.png'
    # plt.savefig(out_fn,bbox_inches='tight')
    # plt.close()

def compute_support(ref_genes, query_exons): #compute the support of all transcripts in ref_genes
    # get overlapping genes:
    # query_exons=junctions_from_read(read)
    query_len = sum([e[1]-e[0] for e in query_exons])
    query_nsj = len(query_exons)*2-2
    support = {k: list() for k in ['ref_gene', 'ref_transcript',
                                   'ref_tss', 'ref_pas', 'ref_len', 'ref_nSJ', 'exI', 'sjI']}
    # check strand??
    for gene in ref_genes:
        for t_id, tr in gene.data['transcripts'].items():
            db_exons=tr['exons']
            support['ref_gene'].append(gene.data['Name'])
            support['ref_transcript'].append(t_id)
            support['ref_len'].append(sum([e[1]-e[0] for e in db_exons]))
            support['ref_nSJ'].append(len(db_exons)*2-2)
            support['ref_tss'].append(
                db_exons[0][0] if gene.data['strand'] == '+' else db_exons[-1][1])
            support['ref_pas'].append(
                db_exons[0][0] if gene.data['strand'] == '-' else db_exons[-1][1])
            sji, exi = get_intersects(db_exons, query_exons)
            support["sjI"].append(sji)
            support['exI'].append(exi)
    support['exIoU'] = [i/(query_len+db-i)
                        for i, db in zip(support['exI'],  support['ref_len'])]
    support['sjIoU'] = [i/(query_nsj+db-i) if (query_nsj+db-i) >
                        0 else 1 for i, db in zip(support['sjI'],  support['ref_nSJ'])]
    return support

def get_splice_type(ref, alt, is_reversed=False):
    if len(ref) == 0:
        return(['novel/unknown'])
    types = ['alternative_donor', 'alternative_acceptor', 'alternative_promoter', 'alternative_polyA',
             'truncated5', 'truncated3', 'exon_skipping', 'novel_exon', 'gapped_exon', 'retained_intron']
    types = {t: [] for t in types}
    relation = get_relation(ref, alt)
    # alt exons that do not overlap any ref
    first = next(i for i, rel in enumerate(relation) if rel) #first ref exon that overlaps an alt exon
    last=len(relation)-next(iter([i for i,r in enumerate(reversed(relation)) if r]))-1
    novel = set(range(len(alt)))-{r[0] for x in relation for r in x}
    if 0 in novel:
        types['alternative_polyA' if is_reversed else 'alternative_promoter'].append(
            '{}-{}'.format(*alt[0]))
    if len(alt)-1 in novel:
        types['alternative_promoter' if is_reversed else 'alternative_polyA'].append(
            '{}-{}'.format(*alt[-1]))
    if len(novel - {0, len(alt)-1}) > 0:
        # todo: distinguish completely novel from novel for this ref_transcript
        types['novel_exon'] += ['{}-{}'.format(*e)
                                for i, e in enumerate(alt) if i in novel-{0, len(alt)-1}]
    # if 0 not in novel and not types['novel_exon'] and len(relation[0]) == 0:

    all_match = (len(novel) == 0 and 
        all([len(rel) == 1 for rel in relation[first:(last+1)]]) and 
        all([rel[0][1] == 3 for rel in relation[(first+1):last]]) and 
        (relation[first][0][1] & 1) and 
        (relation[last][0][1] & 2))
    if first > 0 and all_match:
        types['truncated3' if is_reversed else 'truncated5'].append(
            '{}-{}'.format(*alt[0]))
    # if len(alt)-1 not in novel and not types['novel_exon'] and len(relation[-1]) == 0:
    if last < len(relation)-1 and all_match:
        types['truncated5' if is_reversed else 'truncated3'].append(
            '{}-{}'.format(*alt[-1]))
    for i in range(first, last+1):
        rel = relation[i]
        if len(rel) > 1:  # more than one alt exon for a ref exon
            types['gapped_exon'].append(
                "~".join(['{}-{}'.format(*e) for e in ref]))
        if rel and i > first and relation[i-1] and rel[0][0] == relation[i-1][-1][0]:
            types['retained_intron'].append('{}-{}'.format(*alt[rel[0][0]]))
        # fst splice site does not correspond to any alt exon
        if rel and i < last and rel[-1][1] & 1 == 0 and i < last:
            delta = alt[rel[-1][0]][1]-ref[i][1]
            types['alternative_donor' if is_reversed else 'alternative_acceptor'].append(
                (alt[rel[-1][0]][1], delta))
        if rel and rel[0][1] & 2 == 0 and i > first:
            delta = alt[rel[0][0]][0]-ref[i][0]
            types['alternative_acceptor' if is_reversed else 'alternative_donor'].append(
                (alt[rel[0][0]][0], delta))
        if not rel and i > first and i < last:  # exon skipping
            pre=next(((j,relation[j][-1]) for j in reversed(range(i)) if relation[j]),[]) #first ref exon that overlaps an alt exon
            suc=next(((j,relation[j][0]) for j in range(i+1, len(relation)) if relation[j]),[]) #first ref exon that overlaps an alt exon
            # only predecessing as well as successing splice site is identical
            #if pre and suc and pre[1][1] & 1 and suc[1][1] & 2:
            types['exon_skipping'].append(
                    '{}~{}'.format(alt[pre[1][0]][1], alt[suc[1][0]][0]))
    # return only types that are present
    return {k: v for k, v in types.items() if v}

def get_relation(exons1, exons2):
    # returns a list of length of exons1.
    # each element represents a exon from set 1 and contains a list of pairs, denoting corresponding exons from set 2
    # the first element of the pair is the index of the corresponding set 2 exon,
    # the second element is in [0,1,2,3] and encodes the splice correspondence: 0 overlap only, 1: same splice donor, 2: same splice acceptor, 3: identical
    relation = [[] for _ in exons1]
    enum2 = enumerate(exons2)
    i, e2 = next(enum2)
    for j, e1 in enumerate(exons1):
        while e1[1] > e2[0]:
            ol = overlap(e1, e2)-1
            if ol >= 0:
                if e1[1] == e2[1]:
                    ol += 1
                if e1[0] == e2[0]:
                    ol += 2
                relation[j].append((i, ol))
            if e2[1]<=e1[1]:
                try:
                    i, e2 = next(enum2)
                except StopIteration:  # exons2 is at end
                    return relation
            else: break
    return relation

def overlap(r1, r2):
    # assuming start < end
    if r1[1] < r2[0] or r2[1] < r1[0]:
        return False
    else:
        return True

def get_intersects(tr1, tr2):
    tr1_enum = enumerate(tr1)
    try:
        j,tr1_exon = next(tr1_enum)
    except StopIteration:
        return 0, 0
    sjintersect = 0
    intersect = 0
    for i, tr2_exon in enumerate(tr2):
        while tr1_exon[0] < tr2_exon[1]:
            if tr2_exon[0] == tr1_exon[0] and i > 0 and j>0:  # neglegt TSS and polyA
                sjintersect += 1
            if tr2_exon[1] == tr1_exon[1] and i < len(tr2)-1 and j <len(tr1)-1:
                sjintersect += 1
            if overlap(tr1_exon, tr2_exon):
                # region intersect
                i_end = min(tr1_exon[1], tr2_exon[1])
                i_start = max(tr1_exon[0], tr2_exon[0])
                intersect += (i_end-i_start)
            try:
                j,tr1_exon = next(tr1_enum)
            except StopIteration:  # tr1 is at end
                return sjintersect, intersect
    # tr2 is at end
    return sjintersect, intersect

def junctions_from_read(read):
    exons = list([[read.reference_start, read.reference_start]])
    for cigar in read.cigartuples:
        if cigar[0] == 3:  # N ->  Splice junction
            pos = exons[-1][1]+cigar[1]
            exons.append([pos, pos])
        elif cigar[0] in (0, 2, 7, 8):  # MD=X -> move forward on reference
            exons[-1][1] += cigar[1]
    return exons

def collapse_isoforms(bam_fn='/project/42/pacbio/hecatos_isoseq/05-gmap/all_merged_hq_isoforms.bam'):
    genes = dict()
    gene_number = 0
    isoform_number = 0
    with AlignmentFile(bam_fn, "rb") as align:
        stats = align.get_index_statistics()
        # try catch if sam/ no index
        total_reads = sum([s.mapped for s in stats])
        for i, read in tqdm(enumerate(align.fetch()), total=total_reads):
            tr_name = read.query_name[:read.query_name.find('|')]
            chr = read.reference_name
            strand = '-' if read.is_reverse else '+'
            genes.setdefault(chr, IntervalTree())
            exons = junctions_from_read(read)
            pos = exons[0][0], exons[-1][1]
            for gene in genes[chr].overlap(*pos):
                same_gene = False
                if gene.data[0] == strand:
                    for j, tr in enumerate(gene.data[1]):
                        if splice_identical(tr, exons):
                            tr[0][0] = min(pos[0], tr[0][0])
                            tr[-1][1] = max(exons[-1][1], pos[1])
                            gene.data[3][j].append(pos)
                            same_gene = True
                            break
                        # same gene: at least one shared splice site
                        elif any([r[1] > 0 for x in get_relation(exons, tr) for r in x if r]):
                            same_gene = True  # continue looking for a better match
                            
                    else:
                        if same_gene:
                            gene.data[1].append(exons)
                            gene.data[3].append([pos])
                            isoform_number += 1
                if same_gene:
                    break
            else:  # new gene
                genes[chr].addi(
                    *pos, (strand, [exons], 'PB_{}'.format(gene_number), [[pos]]))
                gene_number += 1
                isoform_number += 1

    print('collapsed {} pacbio transcripts to {} genes and {} isoforms'.format(
        total_reads, gene_number, isoform_number))
    return genes

def merge_interval(e1, e2, must_overlap=False):
    if must_overlap and not overlap(e1, e2):
        raise ValueError('intervals {} and {} do not overlap'.format(e1, e2))
    return (min(e1[0], e2[0]), max(e1[1], e2[1]))

def shared_exons(tr1, tr2):
    rel = get_relation(tr1, tr2)
    same = {r[0] for x in rel for r in x if r[1] == 3}

    return len(same)

def splice_identical(tr1, tr2):
    if len(tr1) != len(tr2):  # different number of exons
        return False
    if len(tr1) == 1 and overlap(tr1[0], tr2[0]):  # single exon genes
        return True
    if tr1[0][1] != tr2[0][1] or tr1[-1][0] != tr2[-1][0]:  # check first and last exons
        return False
    for e1, e2 in zip(tr1[1:-1], tr2[1:-1]):  # check other exons
        if e1[0] != e2[0] or e1[1] != e2[1]:
            return False
    return True


if __name__ == '__main__':
    my_gene=Interval(225494814, 225514661, {'strand': '-', 'Name': 'PB_3316', 'source': {'PB_3316_1': ['all_merged_transcript/21239|full_length_coverage=3|length=3588|num_subreads=57']}, 'transcripts': {'PB_3316_1': {'exons': [(225494814, 225497812), (225494814, 225497812), [225494814, 225497812], [225498346, 225498404], [225500991, 225501070], [225504990, 225505053], [225507950, 225508017], [225511810, 225511859], [225512656, 225512714], [225512870, 225513016], [225514595, 225514661]], 'source': ['all_merged_transcript/58366|full_length_coverage=11|length=2276|num_subreads=60', 'all_merged_transcript/21999|full_length_coverage=2|length=3522|num_subreads=16', 'all_merged_transcript/21239|full_length_coverage=3|length=3588|num_subreads=57'], 'support': {'ref_gene': 'ENAH', 'ref_transcript': 'XM_017001748.1', 'ref_tss': 225653878, 'ref_pas': 225492835, 'ref_len': 7627, 'ref_nSJ': 28, 'exI': 3584, 'sjI': 16, 'exIoU': 0.26308448946634366, 'sjIoU': 0.5, 'sType': {'gapped_exon': ['225492835-225497812~225498346-225498404~225500991-225501070~225504990-225505053~225507950-225508017~225511810-225511859~225512656-225512714~225512870-225513016~225514595-225514900~225517195-225517306~225519197-225519565~225530553-225530638~225554905-225555083~225567248-225567414~225652961-225653878']}}}, 'PB_3316_3': {'exons': [(225495476, 225497812), [225495476, 225497812], [225498346, 225498404], [225500991, 225501054], [225504998, 225505051]], 'source': ['all_merged_transcript/55684|full_length_coverage=5|length=2379|num_subreads=60', 'all_merged_transcript/51489|full_length_coverage=4|length=2511|num_subreads=60'], 'support': {'ref_gene': 'ENAH', 'ref_transcript': 'XM_017001752.1', 'ref_tss': 225653878, 'ref_pas': 225492835, 'ref_len': 7453, 'ref_nSJ': 24, 'exI': 2457, 'sjI': 4, 'exIoU': 0.2496443812233286, 'sjIoU': 0.14285714285714285, 'sType': {'alternative_promoter': ['225504998-225505051'], 'gapped_exon': ['225492835-225497812~225498346-225498404~225500991-225501070~225507950-225508017~225511810-225511859~225512656-225512714~225512870-225513016~225514595-225514900~225519197-225519565~225530553-225530638~225554905-225555083~225567248-225567414~225652961-225653878']}}}, 'PB_3316_4': {'exons': [(225495477, 225497812), [225495477, 225497812], [225498346, 225498404], [225500991, 225501070], [225504990, 225505053], [225507950, 225508017], [225511810, 225511859], [225512656, 225512714], [225512870, 225513016], [225514595, 225514900], [225519197, 225519565], [225530553, 225530638], [225554905, 225555083], [225567248, 225567414], [225652685, 225652843]], 'source': ['all_merged_transcript/20134|full_length_coverage=2|length=3684|num_subreads=23', 'all_merged_transcript/13293|full_length_coverage=2|length=4115|num_subreads=60'], 'support': {'ref_gene': 'ENAH', 'ref_transcript': 'NM_001008493.2', 'ref_tss': 225653143, 'ref_pas': 225486828, 'ref_len': 13175, 'ref_nSJ': 28, 'exI': 4115, 'sjI': 26, 'exIoU': 0.26531270148291425, 'sjIoU': 0.8666666666666667, 'sType': {'exon_skipping': ['225514900~225519197'], 'gapped_exon': ['225486828-225497812~225498346-225498404~225500991-225501070~225504990-225505053~225507950-225508017~225511810-225511859~225512656-225512714~225512870-225513016~225514595-225514900~225517195-225517306~225519197-225519565~225530553-225530638~225554905-225555083~225567248-225567414~225652685-225653143']}}}, 'PB_3316_6': {'exons': [(225495477, 225497812), [225495477, 225497812], [225498346, 225498404], [225500991, 225501070], [225507950, 225508017], [225511810, 225511859], [225512656, 225512714], [225512870, 225513016], [225514595, 225514900], [225517195, 225517306], [225519197, 225519565], [225530553, 225530579]], 'source': ['all_merged_transcript/24488|full_length_coverage=2|length=3448|num_subreads=9', 'all_merged_transcript/20508|full_length_coverage=2|length=3602|num_subreads=31'], 'support': {'ref_gene': 'ENAH', 'ref_transcript': 'NM_018212.5', 'ref_tss': 225653143, 'ref_pas': 225486828, 'ref_len': 13112, 'ref_nSJ': 26, 'exI': 3602, 'sjI': 20, 'exIoU': 0.23318443710752898, 'sjIoU': 0.7142857142857143, 'sType': {'gapped_exon': ['225486828-225497812~225498346-225498404~225500991-225501070~225507950-225508017~225511810-225511859~225512656-225512714~225512870-225513016~225514595-225514900~225517195-225517306~225519197-225519565~225530553-225530638~225554905-225555083~225567248-225567414~225652685-225653143']}}}, 'PB_3316_7': {'exons': [(225495477, 225497812), [225495477, 225497812], [225498346, 225498404], [225500991, 225501070], [225504990, 225505053], [225507950, 225508017], [225511810, 225511859], [225512656, 225512714], [225512870, 225513016], [225514595, 225514900], [225519197, 225519567]], 'source': ['all_merged_transcript/33413|full_length_coverage=2|length=3049|num_subreads=60', 'all_merged_transcript/23227|full_length_coverage=4|length=3559|num_subreads=36'], 'support': {'ref_gene': 'ENAH', 'ref_transcript': 'XM_017001748.1', 'ref_tss': 225653878, 'ref_pas': 225492835, 'ref_len': 7627, 'ref_nSJ': 28, 'exI': 3528, 'sjI': 18, 'exIoU': 0.35407466880770777, 'sjIoU': 0.6, 'sType': {'exon_skipping': ['225514900~225519197'], 'gapped_exon': ['225492835-225497812~225498346-225498404~225500991-225501070~225504990-225505053~225507950-225508017~225511810-225511859~225512656-225512714~225512870-225513016~225514595-225514900~225517195-225517306~225519197-225519565~225530553-225530638~225554905-225555083~225567248-225567414~225652961-225653878']}}}, 'PB_3316_10': {'exons': [(225495477, 225497812), [225495477, 225497812], [225498346, 225498404], [225500991, 225501070], [225504990, 225505053], [225507950, 225508017], [225511810, 225511859], [225512656, 225512714], [225512870, 225513016], [225514595, 225514615]], 'source': ['all_merged_transcript/43956|full_length_coverage=3|length=2694|num_subreads=49', 'all_merged_transcript/38695|full_length_coverage=4|length=2876|num_subreads=52'], 'support': {'ref_gene': 'ENAH', 'ref_transcript': 'XM_017001748.1', 'ref_tss': 225653878, 'ref_pas': 225492835, 'ref_len': 7627, 'ref_nSJ': 28, 'exI': 2875, 'sjI': 16, 'exIoU': 0.28859666733587636, 'sjIoU': 0.5333333333333333, 'sType': {'gapped_exon': ['225492835-225497812~225498346-225498404~225500991-225501070~225504990-225505053~225507950-225508017~225511810-225511859~225512656-225512714~225512870-225513016~225514595-225514900~225517195-225517306~225519197-225519565~225530553-225530638~225554905-225555083~225567248-225567414~225652961-225653878']}}}, 'PB_3316_12': {'exons': [(225495477, 225497812), [225495477, 225497812], [225498346, 225498404], [225500991, 225501070], [225504990, 225505053], [225507950, 225507976]], 'source': ['all_merged_transcript/46746|full_length_coverage=4|length=2647|num_subreads=60', 'all_merged_transcript/48734|full_length_coverage=4|length=2562|num_subreads=60'], 'support': {'ref_gene': 'ENAH', 'ref_transcript': 'XM_017001748.1', 'ref_tss': 225653878, 'ref_pas': 225492835, 'ref_len': 7627, 'ref_nSJ': 28, 'exI': 2561, 'sjI': 8, 'exIoU': 0.2570768921903232, 'sjIoU': 0.26666666666666666, 'sType': {'gapped_exon': ['225492835-225497812~225498346-225498404~225500991-225501070~225504990-225505053~225507950-225508017~225511810-225511859~225512656-225512714~225512870-225513016~225514595-225514900~225517195-225517306~225519197-225519565~225530553-225530638~225554905-225555083~225567248-225567414~225652961-225653878']}}}, 'PB_3316_16': {'exons': [(225495479, 225497812), [225495479, 225497812], [225498346, 225498404], [225500991, 225501070], [225507950, 225508017], [225511810, 225511859], [225512656, 225512714], [225512870, 225513016], [225514595, 225514900], [225519197, 225519565], [225530553, 225530638], [225554905, 225555011]], 'source': ['all_merged_transcript/29643|full_length_coverage=4|length=3226|num_subreads=60', 'all_merged_transcript/19330|full_length_coverage=2|length=3657|num_subreads=42'], 'support': {'ref_gene': 'ENAH', 'ref_transcript': 'XM_017001752.1', 'ref_tss': 225653878, 'ref_pas': 225492835, 'ref_len': 7453, 'ref_nSJ': 24, 'exI': 3654, 'sjI': 20, 'exIoU': 0.37339055793991416, 'sjIoU': 0.7692307692307693, 'sType': {'gapped_exon': ['225492835-225497812~225498346-225498404~225500991-225501070~225507950-225508017~225511810-225511859~225512656-225512714~225512870-225513016~225514595-225514900~225519197-225519565~225530553-225530638~225554905-225555083~225567248-225567414~225652961-225653878']}}}, 'PB_3316_18': {'exons': [[225495552, 225497812], [225498346, 225498404], [225500991, 225501070], [225507950, 225507992]], 'source': ['all_merged_transcript/52570|full_length_coverage=2|length=2441|num_subreads=36'], 'support': {'ref_gene': 'ENAH', 'ref_transcript': 'XM_017001752.1', 'ref_tss': 225653878, 'ref_pas': 225492835, 'ref_len': 7453, 'ref_nSJ': 24, 'exI': 2439, 'sjI': 6, 'exIoU': 0.32725077150140885, 'sjIoU': 0.25, 'sType': {'truncated5': ['225507950-225507992']}}}, 'PB_3316_19': {'exons': [[225495673, 225497812], [225498346, 225498404], [225500991, 225501070], [225504990, 225505053], [225507950, 225507978]], 'source': ['all_merged_transcript/55818|full_length_coverage=2|length=2368|num_subreads=36'], 'support': {'ref_gene': 'ENAH', 'ref_transcript': 'XM_017001748.1', 'ref_tss': 225653878, 'ref_pas': 225492835, 'ref_len': 7627, 'ref_nSJ': 28, 'exI': 2367, 'sjI': 8, 'exIoU': 0.3103448275862069, 'sjIoU': 0.2857142857142857, 'sType': {'truncated5': ['225507950-225507978']}}}}, 'splice_graph': IntervalTree([Interval(225494814, 225495476, (set(), {225495476})), Interval(225495476, 225495477, ({225495476}, {225495477})), Interval(225495477, 225495479, ({225495477}, {225495479})), Interval(225495479, 225495552, ({225495479}, {225495552})), Interval(225495552, 225495673, ({225495552}, {225495673})), Interval(225495673, 225497812, ({225495673}, {225498346})), Interval(225498346, 225498404, ({225497812}, {225500991})), Interval(225500991, 225501054, ({225498404}, {225504998, 225501070, 225501054})), Interval(225501054, 225501070, ({225501054}, {225504990, 225507950})), Interval(225504990, 225504998, ({225501070}, {225504998})), Interval(225504998, 225505051, ({225504998, 225501054}, {225505051, 225505053})), Interval(225505051, 225505053, ({225505051}, {225507950})), Interval(225507950, 225507976, ({225505053, 225501070}, {225507976, 225508017})), Interval(225507976, 225507978, ({225507976}, {225507992})), Interval(225507978, 225507992, (set(), {225508017})), Interval(225507992, 225508017, (set(), {225511810})), Interval(225511810, 225511859, ({225508017}, {225512656})), Interval(225512656, 225512714, ({225511859}, {225512870})), Interval(225512870, 225513016, ({225512714}, {225514595})), Interval(225514595, 225514615, ({225513016}, {225514661, 225514615})), Interval(225514615, 225514661, ({225514615}, {225514661})), Interval(225514661, 225514900, ({225514661}, {225517195, 225519197})), Interval(225517195, 225517306, ({225514900}, {225519197})), Interval(225519197, 225519565, ({225517306, 225514900}, {225530553, 225519565})), Interval(225519565, 225519567, ({225519565}, set())), Interval(225530553, 225530579, ({225519565}, {225530579, 225530638})), Interval(225530579, 225530638, ({225530579}, {225554905})), Interval(225554905, 225555011, ({225530638}, {225555083})), Interval(225555011, 225555083, (set(), {225567248})), Interval(225567248, 225567414, ({225555083}, {225652685})), Interval(225652685, 225652843, ({225567414}, set()))])})
    collapse_transcript_of_gene(my_gene, 5,1000,100, False)
    [i[1]-i[0] for i in sorted(my_gene.data['splice_graph'])]

    arg_parser = argparse.ArgumentParser(
        description='compare isoseq alignment to refseq')
    arg_parser.add_argument('-i', '--bam_fn',    help='input bam file', type=str, default='/project/42/pacbio/hecatos_isoseq/05-gmap/all_merged_hq_isoforms.bam')
    arg_parser.add_argument('-r', '--refseq_fn', help='refseq gtf file', type=str)
    arg_parser.add_argument('-o', '--outsuffix', help='suffix for output file', type=str, default='_refseq_tab.txt')
    arg_parser.add_argument('-a', '--analyse',   help='analyse table', type=str)
    args=arg_parser.parse_args()
    #isoseq=import_bam_transcripts(args.bam_fn, n_lines=1000)
    #collapse_transcripts(isoseq)

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