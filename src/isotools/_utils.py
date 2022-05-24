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
                    pbar.update(1)
                    continue
                readl = len(read.query_qualities)
                if len_bins is not None:
                    len_i = next((i for i, th in enumerate(len_bins) if readl < th), len(len_bins))
                error_rate = sum([10**(-q/10) for q in read.query_qualities])/readl*100
                q_i = next((i for i, th in enumerate(qual_bins) if error_rate < th), len(qual_bins))
                qual[q_i, len_i] += 1
                pbar.update(1)
                i += 1
                if i >= n:
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
    if len(tr1) == 1 and has_overlap(tr1[0], tr2[0]):  # single exon genes
        return True
    if tr1[0][1] != tr2[0][1] or tr1[-1][0] != tr2[-1][0]:  # check first and last exons
        return False
    for e1, e2 in zip(tr1[1:-1], tr2[1:-1]):  # check other exons
        if e1[0] != e2[0] or e1[1] != e2[1]:
            return False
    return True


def has_overlap(r1, r2):
    "check the overlap of two intervals"
    # assuming start < end
    if r1[1] < r2[0] or r2[1] < r1[0]:
        return False
    else:
        return True


def get_overlap(r1, r2):
    "check the overlap of two intervals"
    # assuming start < end
    return max(0, min(r1[1], r2[1]) - max(r1[0], r2[0]))


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
            if has_overlap(tr1_exon, tr2_exon):
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


def _interval_dist(a, b):
    '''compute the distance between two intervals a and b.'''
    return max([a[0], b[0]])-min([a[1], b[1]])


def _filter_event(coverage, event, min_total=100, min_alt_fraction=.1):
    '''
    return True if the event satisfies the filter conditions and False otherwise

    :param coverage: 1D array of counts per transcript
    :param event: Event obtained from .find_splice_bubbles()
    :param min_total: The minimum total number of reads for an event to pass the filter
    :type min_total: int
    :param min_alt_fraction: The minimum fraction of read supporting the alternative
    :type min_alt_frction: float'''

    tr_IDs = event[0]+event[1]
    tot_cov = coverage[tr_IDs].sum()

    if tot_cov < min_total:
        return False

    pri_cov = coverage[event[0]].sum()
    alt_cov = coverage[event[1]].sum()
    frac = min(pri_cov, alt_cov)/tot_cov

    if frac < min_alt_fraction:
        return False

    return True


def _get_exonic_region(transcripts):
    e_starts = iter(sorted([e[0] for tr in transcripts for e in tr['exons']]))
    e_ends = iter(sorted([e[1] for tr in transcripts for e in tr['exons']]))
    exon_reg = [[next(e_starts), next(e_ends)]]
    for next_start in e_starts:
        if next_start <= exon_reg[-1][1]:
            exon_reg[-1][1] = next(e_ends)
        else:
            exon_reg.append([next_start, next(e_ends)])
    return exon_reg


def _get_overlap(exons, transcripts):
    '''Compute the exonic overlap of a new transcript with the segment graph.
    Avoids the computation of segment graph, which provides the same functionality.

    :param exons: A list of exon tuples representing the new transcript
    :type exons: list
    :return: boolean array indicating whether the splice site is contained or not'''
    if not transcripts:
        return 0
    # 1) get exononic regions in transcripts
    exon_reg = _get_exonic_region(transcripts)
    # 2) find overlap of exonic regions with exons
    ol = 0
    i = 0
    for e in exons:
        while(exon_reg[i][1] < e[0]):  # no overlap, go on
            i += 1
            if i == len(exon_reg):
                return ol
        while(exon_reg[i][0] < e[1]):
            i_end = min(e[1], exon_reg[i][1])
            i_start = max(e[0], exon_reg[i][0])
            ol += (i_end - i_start)
            if exon_reg[i][1] > e[1]:  # might overlap with next exon
                break
            i += 1
            if i == len(exon_reg):
                return ol
    return ol


def _find_splice_sites(sj, transcripts):
    '''Checks whether the splice sites of a new transcript are present in the set of transcripts.
    Avoids the computation of segment graph, which provides the same functionality.

    :param sj: A list of 2 tuples with the splice site positions
    :type exons: list
    :param transcripts: transcripts to scan
    :return: boolean array indicating whether the splice site is contained or not'''

    sites = np.zeros((len(sj)) * 2, dtype=bool)
    # check exon ends
    splice_junction_starts = {}
    splice_junction_ends = {}
    for i, ss in enumerate(sj):
        splice_junction_starts.setdefault(ss[0], []).append(i)
        splice_junction_ends.setdefault(ss[1], []).append(i)

    tr_list = [iter(tr['exons'][:-1]) for tr in transcripts if len(tr['exons']) > 1]
    current = [next(tr) for tr in tr_list]
    for sjs, idx in sorted(splice_junction_starts.items()):  # splice junction starts, sorted by position
        for j, tr_iter in enumerate(tr_list):
            try:
                while(sjs > current[j][1]):
                    current[j] = next(tr_iter)
                if current[j][1] == sjs:
                    for i in idx:
                        sites[i * 2] = True
                    break
            except StopIteration:
                continue
    # check exon starts
    tr_list = [iter(tr['exons'][1:]) for tr in transcripts if len(tr['exons']) > 1]
    current = [next(tr) for tr in tr_list]
    for sje, idx in sorted(splice_junction_ends.items()):  # splice junction ends, sorted by position
        for j, tr_iter in enumerate(tr_list):
            try:
                while(sje > current[j][0]):
                    current[j] = next(tr_iter)
                if current[j][0] == sje:
                    for i in idx:
                        sites[i * 2 + 1] = True
                    break
            except StopIteration:
                continue
    return sites


def _corrected_log2OR(con_tab):
    con_tab_copy = np.zeros((2, 2), dtype=float)

    for m, n in itertools.product(range(2), range(2)):
        if con_tab[n, m] == 0:
            con_tab_copy[n, m] = 10**-9
        else:
            con_tab_copy[n, m] = con_tab[n, m]
    log2OR = np.log2((con_tab_copy[0, 0]*con_tab_copy[1, 1])) - np.log2((con_tab_copy[0, 1]*con_tab_copy[1, 0]))
    return log2OR
