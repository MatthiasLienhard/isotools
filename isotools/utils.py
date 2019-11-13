
import warnings
from itertools import combinations
import itertools
import numpy as np
import matplotlib.pyplot as plt


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

def junctions_from_read(read):
    exons = list([[read.reference_start, read.reference_start]])
    for cigar in read.cigartuples:
        if cigar[0] == 3:  # N ->  Splice junction
            pos = exons[-1][1]+cigar[1]
            exons.append([pos, pos])
        elif cigar[0] in (0, 2, 7, 8):  # MD=X -> move forward on reference
            exons[-1][1] += cigar[1]
    return exons

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


def overlap(r1, r2):
    # assuming start < end
    if r1[1] < r2[0] or r2[1] < r1[0]:
        return False
    else:
        return True


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
