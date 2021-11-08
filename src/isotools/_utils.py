from pysam import AlignmentFile
import numpy as np
import pandas as pd
import itertools
import re
from tqdm import tqdm
import builtins
import logging
logger = logging.getLogger('isotools')

cigar = 'MIDNSHP=XB'
cigar_lup = {c: i for i, c in enumerate(cigar)}

compl = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def rc(seq):
    '''reverse complement of seq'''
    return ''.join(reversed([compl[c] if c in compl else 'N' for c in seq]))


def get_error_rate(bam_fn, n=1000):
    qual = 0
    total_len = 0
    with AlignmentFile(bam_fn, "rb", check_sq=False) as align:
        if n is None:
            stats = align.get_index_statistics()
            n = sum([s.mapped for s in stats])
        with tqdm(total=n, unit=' reads') as pbar:
            for i, read in enumerate(align):
                total_len += len(read.query_qualities)
                qual += sum([10**(-q/10) for q in read.query_qualities])
                pbar.update(1)
                if i+1 >= n:
                    break
    return (qual/total_len)*100


def basequal_hist(bam_fn, qual_bins=10**(np.linspace(-7, 0, 30)), len_bins=None, n=10000):
    '''compute summary base quality statistics from base file, optionally dependent on read length'''
    n_len_bins = 1 if len_bins is None else len(len_bins)+1
    qual = np.zeros((len(qual_bins)+1, n_len_bins), dtype=int)
    len_i = 0
    i = 0
    with AlignmentFile(bam_fn, "rb") as align:
        if n is None:
            stats = align.get_index_statistics()
            n = sum([s.mapped for s in stats])
        with tqdm(total=n, unit=' reads') as pbar:
            for read in align:
                if read.query_qualities is None:
                    continue
                readl = len(read.query_qualities)
                if len_bins is not None:
                    len_i = next((i for i, th in enumerate(len_bins) if readl < th), len(len_bins))
                error_rate = sum([10**(-q/10) for q in read.query_qualities])/readl*100
                q_i = next((i for i, th in enumerate(qual_bins) if error_rate < th), len(qual_bins))
                qual[q_i, len_i] += 1
                pbar.update(1)
                i += 1
                if i+1 >= n:
                    break
    idx = [f'<{th:.2E} %' for th in qual_bins]+[f'>={qual_bins[-1]:.2E} %']
    if len_bins is None:
        return pd.Series(qual[:, 0], index=idx)
    col = [f'<{th/1000:.1f} kb' for th in len_bins]+[f'>={len_bins[-1]/1000:.1f} kb']
    return pd.DataFrame(qual, index=idx, columns=col)


def pairwise(iterable):  # e.g. usefull for enumerating introns
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def cigar_string2tuples(cigarstring):
    res = re.findall(f'(\\d+)([{cigar}]+)', cigarstring)
    return tuple((cigar_lup[c], int(n)) for n, c in res)


def junctions_from_cigar(cigartuples, offset):
    'returns the exon positions'
    exons = list([[offset, offset]])
    for cigar in cigartuples:
        if cigar[0] == 3:  # N ->  Splice junction
            pos = exons[-1][1]+cigar[1]
            if exons[-1][0] == exons[-1][1]:
                # delete zero length exons
                # (may occur if insertion within intron, e.g. 10M100N10I100N10M)
                del exons[-1]
            exons.append([pos, pos])
        elif cigar[0] in (0, 2, 7, 8):  # MD=X -> move forward on reference
            exons[-1][1] += cigar[1]
    if exons[-1][0] == exons[-1][1]:  # delete 0 length exons at the end
        del exons[-1]
    return exons


def is_same_gene(tr1, tr2, spj_iou_th=0, reg_iou_th=.5):
    'checks whether tr1 and tr2 are the same gene by calculating intersection over union of the intersects'
    # current default definition of "same gene": at least one shared splice site
    # or more than 50% exonic overlap
    spj_i, reg_i = get_intersects(tr1, tr2)
    total_spj = (len(tr1)+len(tr2)-2)*2
    spj_iou = spj_i/(total_spj-spj_i) if total_spj > 0 else 0
    if spj_iou > spj_iou_th:
        return True
    total_len = sum([e[1]-e[0] for e in tr2+tr1])
    reg_iou = reg_i/(total_len-reg_i)
    if reg_iou > reg_iou_th:
        return True
    return False


def splice_identical(tr1, tr2):
    # all splice sites are equal
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


def overlap(r1, r2):
    "check the overlap of two intervals"
    # assuming start < end
    if r1[1] < r2[0] or r2[1] < r1[0]:
        return False
    else:
        return True


def get_intersects(tr1, tr2):
    "get the number of intersecting splice sites and intersecting bases of two transcripts"
    tr1_enum = enumerate(tr1)
    try:
        j, tr1_exon = next(tr1_enum)
    except StopIteration:
        return 0, 0
    sjintersect = 0
    intersect = 0
    for i, tr2_exon in enumerate(tr2):
        while tr1_exon[0] < tr2_exon[1]:
            if tr2_exon[0] == tr1_exon[0] and i > 0 and j > 0:  # neglegt TSS and polyA
                sjintersect += 1
            if tr2_exon[1] == tr1_exon[1] and i < len(tr2)-1 and j < len(tr1)-1:
                sjintersect += 1
            if overlap(tr1_exon, tr2_exon):
                # region intersect
                i_end = min(tr1_exon[1], tr2_exon[1])
                i_start = max(tr1_exon[0], tr2_exon[0])
                intersect += (i_end-i_start)
            try:
                j, tr1_exon = next(tr1_enum)
            except StopIteration:  # tr1 is at end
                return sjintersect, intersect
    # tr2 is at end
    return sjintersect, intersect


def _filter_function(expression):
    'converts a string e.g. "all(x[0]/x[1]>3) " into a function'
    # extract argument names
    f = eval(f'lambda: {expression}')
    args = [n for n in f.__code__.co_names if n not in dir(builtins)]

    # potential issue: g.coverage gets detected as ["g", "coverage"], e.g. coverage is added. Probably not causing trubble
    return eval(f'lambda {",".join([arg+"=None" for arg in args]+["**kwargs"])}: bool({expression})\n', {}, {}), args
