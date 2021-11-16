import numpy as np
import logging
from sortedcontainers import SortedDict  # for SpliceGraph
from ._utils import pairwise, overlap
from .decorators import deprecated, experimental
from typing import Union

logger = logging.getLogger('isotools')


class SegmentGraph():
    '''Segment Graph Implementation

    Nodes in the Segment Graph represent disjoint exonic bins (aka segments) and have start (genomic 5'), end (genomic 3'),
    and a dict of successors and predesessors (edges)
    Edges link two exonic bins x and y that follow successively in a transcript.
    They represent either introns (if x.end<y.start) or connect exonic bins of the same exon (x.end==y.start).

    :param transcripts: A list of transcripts, which are lists of exons, which in turn are (start,end) tuples
    :type transcripts: list
    :param strand: the strand of the gene, either "+" or "-"'''

    def __init__(self, transcripts, strand):
        self.strand = strand
        assert strand in '+-', 'strand must be either "+" or "-"'
        open_exons = dict()
        for tr in transcripts:
            for e in tr:
                open_exons[e[0]] = open_exons.get(e[0], 0) + 1
                open_exons[e[1]] = open_exons.get(e[1], 0) - 1
        # sort by value
        boundaries = sorted(open_exons)
        open_count = 0
        self._graph = list()
        for i, start in enumerate(boundaries[:-1]):
            open_count += open_exons[start]
            if open_count > 0:
                self._graph.append(SegGraphNode(start, boundaries[i + 1]))

        # get the links
        start_idx = {node.start: i for i, node in enumerate(self._graph)}
        end_idx = {node.end: i for i, node in enumerate(self._graph)}
        self._tss = [start_idx[e[0][0]] for e in transcripts]  # todo: this is missleading: on - strand this would not be the tss
        self._pas = [end_idx[e[-1][1]] for e in transcripts]

        for i, tr in enumerate(transcripts):
            for j, e in enumerate(tr):
                if start_idx[e[0]] < end_idx[e[1]]:  # psydojuctions within exon
                    self._graph[start_idx[e[0]]].suc[i] = start_idx[e[0]] + 1
                    self._graph[end_idx[e[1]]].pre[i] = end_idx[e[1]] - 1
                    for node_idx in range(start_idx[e[0]] + 1, end_idx[e[1]]):
                        self._graph[node_idx].suc[i] = node_idx + 1
                        self._graph[node_idx].pre[i] = node_idx - 1
                if j < len(tr) - 1:  # real junctions
                    e2 = tr[j + 1]
                    self._graph[end_idx[e[1]]].suc[i] = start_idx[e2[0]]
                    self._graph[start_idx[e2[0]]].pre[i] = end_idx[e[1]]

    def _restore(self, i: int) -> list:  # mainly for testing
        ''' Restore the i_{th} transcript from the Segment graph by traversing from 5' to 3'

        :param i: The index of the transcript to restore
        :type i: int
        :return: A list of exon tuples representing the transcript
        :rtype: list'''
        idx = self._tss[i]
        exons = [[self._graph[idx].start, self._graph[idx].end]]
        while True:
            if idx == self._pas[i]:
                break
            idx = self._graph[idx].suc[i]
            if self._graph[idx].start == exons[-1][1]:  # extend
                exons[-1][1] = self._graph[idx].end
            else:
                exons.append([self._graph[idx].start, self._graph[idx].end])
            # print(exons)

        return exons

    def _restore_reverse(self, i: int) -> list:  # mainly for testing
        ''' Restore the ith transcript from the Segment graph by traversing from 3' to 5'

        :param i: The index of the transcript to restore
        :type i: int
        :return: A list of exon tuples representing the transcript
        :rtype: list'''
        idx = self._pas[i]
        exons = [[self._graph[idx].start, self._graph[idx].end]]
        while True:
            if idx == self._tss[i]:
                break
            idx = self._graph[idx].pre[i]
            if self._graph[idx].end == exons[-1][0]:  # extend
                exons[-1][0] = self._graph[idx].start
            else:
                exons.append([self._graph[idx].start, self._graph[idx].end])
            # print(exons)

        exons.reverse()
        return exons

    def search_transcript(self, exons):
        '''Tests if a transcript (provided as list of exons) is contained in self and return the corresponding transcript indices.

        :param exons: A list of exon tuples representing the transcript
        :type exons: list
        :return: a list of supporting transcript indices
        :rtype: list'''

        # fst special case: exons extends segment graph
        if exons[0][1] <= self[0].start or exons[-1][0] >= self[-1].end:
            return []
        # snd special case: single exon transcript: return all overlapping single exon transcripts form sg
        if len(exons) == 1:
            return [trid for trid, (j1, j2) in enumerate(zip(self._tss, self._pas))
                    if self._is_same_exon(trid, j1, j2) and self[j1].start <= exons[0][1] and self[j2].end >= exons[0][0]]
        # all junctions must be contained and no additional
        tr = set(range(len(self._tss)))
        j = 0
        for i, e in enumerate(exons[:-1]):
            while j < len(self) and self[j].end < e[1]:  # check exon (no junction allowed)
                tr -= set(trid for trid, j2 in self[j].suc.items() if self[j].end != self[j2].start)
                j += 1
            if self[j].end != e[1]:
                return []
            # check junction (must be present)
            tr &= set(trid for trid, j2 in self[j].suc.items() if self[j2].start == exons[i + 1][0])
            j += 1
            if len(tr) == 0:
                return tr

        while j < len(self):  # check last exon (no junction allowed)
            tr -= set(trid for trid, j2 in self[j].suc.items() if self[j].end != self[j2].start)
            j += 1
        return [trid for trid in tr]

    def _is_same_exon(self, tr_nr, j1, j2):
        '''Tests if nodes j1 and j2 belong to same exon in transcript tr_nr.'''
        for j in range(j1, j2):
            if tr_nr not in self[j].suc or self[j].suc[tr_nr] > j + 1 or self[j].end != self[j + 1].start:
                return False
        return True

    def _count_introns(self, tr_nr, j1, j2):
        '''Counts the number of junctions between j1 and j2.'''
        logger.debug('counting introns of transcript %i between nodes %i and %i', tr_nr, j1, j2)
        delta = 0
        if j1 == j2:
            return 0
        assert tr_nr in self[j1].suc, f'transcript {tr_nr} does not contain node {j1}'
        while j1 < j2:
            j_next = self[j1].suc[tr_nr]
            if j_next > j1 + 1 or self[j1].end != self[j1 + 1].start:
                delta += 1
            j1 = j_next
        return delta

    def get_node_matrix(self) -> np.array:
        '''Gets the node matrix representation of the segment graph.'''

        return np.array([[True if tss == j or trid in n.pre else False for j, n in enumerate(self)] for trid, tss in enumerate(self._tss)])

    def find_fragments(self):
        '''Finds all fragments (e.g. transcript contained in other transcripts) in the segment graph.'''
        truncated = set()
        contains = {}
        starts = [set() for _ in range(len(self))]
        for trid, idx in enumerate(self._tss):
            starts[idx].add(trid)
        nodes = self.get_node_matrix()
        for trid, (tss, pas) in enumerate(zip(self._tss, self._pas)):
            if trid in truncated:
                continue
            contains[trid] = {trid2 for trid2, (tss2, pas2) in enumerate(zip(self._tss, self._pas)) if trid2 != trid and tss2 >=
                              tss and pas2 <= pas and all(nodes[trid2, tss2:pas2 + 1] == nodes[trid, tss2:pas2 + 1])}
            truncated.update(contains[trid])  # those are not checked

        fragments = {}
        for big, smallL in contains.items():
            if big not in truncated:
                for trid in smallL:
                    delta1 = self._count_introns(big, self._tss[big], self._tss[trid])
                    delta2 = self._count_introns(big, self._pas[trid], self._pas[big])
                    fragments.setdefault(trid, []).append((big, delta1, delta2) if self.strand == '+' else (big, delta2, delta1))
        return fragments

    def get_alternative_splicing(self, exons, alternative=None):
        '''Compares exons to segment graph and returns list of novel splicing events.

        This function computes the novelty class of the provided transcript compared to (reference annotation) transcripts
        from the segment graph. It returns the "squanti category" (0=FSM,1=ISM,2=NIC,3=NNC,4=Novel gene) and the subcategory.

        :param exons: A list of exon tuples representing the transcript
        :type exons: list
        :return: pair with the squanti category number and the subcategories as list of novel splicing events
            that produce the provided transcript from the transcripts in splce graph
        :rtype: tuple'''

        # returns a tuple
        # the sqanti category: 0=FSM,1=ISM,2=NIC,3=NNC,4=Novel gene
        # subcategories: a list of novel splicing events or splice_identical

        # a list of tuples with (1) gene names and (2) junction numbers covered by other genes (e.g. readthrough fusion)
        if alternative is not None and len(alternative) > 0:
            category = 4
            fusion_exons = {int((i + 1) / 2) for j in alternative for i in j[1]}
            altsplice = {'readthrough fusion': alternative}  # other novel events are only found in the primary reference transcript
        else:
            tr = self.search_transcript(exons)
            if tr:
                return 0, {'FSM': tr}
            category = 1
            altsplice = {}
            fusion_exons = set()

        is_reverse = self.strand == '-'
        j1 = next((j for j, n in enumerate(self) if n.end > exons[0][0]))
        # j1: index of first segment ending after exon start (i.e. first overlapping segment)
        j2 = next((j - 1 for j in range(j1, len(self)) if self[j].start >= exons[0][1]), len(self) - 1)
        # j2: index of last segment starting befor exon end (i.e. last overlapping segment)

        # check truncation at begining (e.g. low position)
        if (len(exons) > 1 and  # no mono exon
                not any(j in self._tss for j in range(j1, j2 + 1)) and  # no tss/pas within exon
                self[j1].start <= exons[0][0]):  # start of first exon is exonic in ref
            j0 = max(self._tss[trid] for trid in self[j1].pre)  # j0 is the closest start node
            if any(self[j].end < self[j + 1].start for j in range(j0, j1)):  # assure there is an intron between closest tss/pas and exon
                end = '3' if is_reverse else '5'
                altsplice.setdefault(f'{end}\' fragment', []).append([self[j0].start, exons[0][0]])  # at start (lower position)

        for i, ex1 in enumerate(exons):
            ex2 = None if i + 1 == len(exons) else exons[i + 1]
            if i not in fusion_exons:  # exon belongs to other gene (read through fusion)
                # finds intron retention (NIC), novel exons, novel splice sites,  novel pas/tss  (NNC)
                exon_altsplice, exon_cat = self._check_exon(j1, j2, i == 0, is_reverse, ex1, ex2)
                category = max(exon_cat, category)
                for k, v in exon_altsplice.items():
                    altsplice.setdefault(k, []).extend(v)
                # find j2: index of last segment starting befor exon2 end (i.e. last overlapping  segment)
                if ex2 is not None:
                    if j2 + 1 < len(self):
                        j1, j2, junction_altsplice = self._check_junction(j1, j2, ex1, ex2)  # finds exon skipping and novel junction (NIC)
                        if junction_altsplice and i + 1 not in fusion_exons:
                            category = max(2, category)
                            for k, v in junction_altsplice.items():
                                altsplice.setdefault(k, []).extend(v)
                    else:
                        j1 = len(self)

        # check truncation at end (e.g. high position)
        if (len(exons) > 1 and
                j2 >= j1 and
                not any(j in self._pas for j in range(j1, j2 + 1)) and  # no tss/pas within exon
                self[j2].end >= exons[-1][1]):  # end of last exon is exonic in ref
            try:
                j3 = min(self._pas[trid] for trid in self[j2].suc)  # j3 is the next end node (pas/tss on fwd/rev)
            except ValueError:
                logger.error('\n'.join([str(exons), str(self._pas), str((j1, j2)), str([(j, n) for j, n in enumerate(self)])]))
                raise
            if any(self[j].end < self[j + 1].start for j in range(j2, j3)):  # assure there is an intron between closest tss/pas and exon
                end = '5' if is_reverse else '3'
                altsplice.setdefault(f'{end}\' fragment', []).append([exons[-1][1], self[j3].end])

        if not altsplice:  # all junctions are contained but not all in one transcript
            altsplice = {'novel combination': []}
            category = 2

        return category, altsplice

    def _check_exon(self, j1, j2, is_first, is_reverse, e, e2=None):
        '''checks whether exon is supported by splice graph between nodes j1 and j2

        :param j1: index of first segment ending after exon start (i.e. first overlapping segment)
        :param j2: index of last segment starting befor exon end (i.e. last overlapping  segment)'''

        logger.debug('exon %s between sg node %s and %s/%s (first=%s,rev=%s,e2=%s)', e, j1, j2, len(self), is_first, is_reverse, e2)
        is_last = e2 is None
        altsplice = {}
        category = 0
        if j1 > j2:  # e is not contained at all  -> novel exon (or TSS/PAS if first/last)
            category = 3
            if is_first or is_last:
                altsplice = {'novel intronic PAS' if is_first == is_reverse else 'novel intronic TSS': [e]}
            else:
                altsplice = {'novel exon': [e]}
            j2 = j1
        elif (is_first and is_last):  # mono-exon
            altsplice['mono-exon'] = []
            category = 1
        else:  # check splice sites
            if self[j1][0] != e[0]:  # first splice site missmatch
                if not is_first:
                    # pos="intronic" if self[j1][0]>e[0] else "exonic"
                    kind = '5' if is_reverse else '3'
                    dist = min((self[j][0] - e[0] for j in range(j1, j2 + 1)), key=abs)  # the distance to next junction
                    altsplice[f"novel {kind}' splice site"] = [(e[0], dist)]
                    category = 3
                elif self[j1][0] > e[0] and not any(j in self._tss for j in range(j1, j2 + 1)):  # exon start is intronic in ref
                    site = 'PAS' if is_reverse else 'TSS'
                    altsplice.setdefault(f'novel exonic {site}', []).append((e[0], self[j1][0]))
                    category = max(2, category)
            if self[j2][1] != e[1]:  # second splice site missmatch
                if not is_last:
                    # pos="intronic" if self[j2][1]<e[1] else "exonic"
                    # TODO: could also be a "novel intron", if the next "first" splice site is also novel.
                    kind = '3' if is_reverse else '5'
                    dist = min((self[j][1] - e[1] for j in range(j1, j2 + 1)), key=abs)  # the distance to next junction
                    altsplice.setdefault(f"novel {kind}' splice site", []).append((e[1], dist))
                    category = 3
                elif self[j2][1] < e[1] and not any(j in self._pas for j in range(j1, j2 + 1)):  # exon end is intronic in ref & not overlapping tss
                    site = 'TSS' if is_reverse else 'PAS'
                    altsplice.setdefault(f'novel exonic {site}', []).append((self[j2][1], e[1]))
                    category = max(2, category)

        # find intron retentions
        if j1 < j2 and any(self[ji + 1].start - self[ji].end > 0 for ji in range(j1, j2)):
            gaps = [ji for ji in range(j1, j2) if self[ji + 1].start - self[ji].end > 0]
            if (gaps
                    and not (is_first and any(j in self._tss for j in range(gaps[-1] + 1, j2)))
                    and not (is_last and any(j in self._pas for j in range(j1, gaps[0] + 1)))):
                ret_introns = []
                troi = set(self[j1].suc.keys()).intersection(self[j2].pre.keys())
                if troi:
                    j = j1
                    while j < j2:
                        nextj = min(js for trid, js in self[j].suc.items() if trid in troi)
                        if self[nextj].start - self[j].end > 0 and any(self[ji + 1].start - self[ji].end > 0 for ji in range(j, nextj)):
                            ret_introns.append((self[j].end, self[nextj].start))
                        j = nextj
                    if ret_introns:
                        altsplice['intron retention'] = ret_introns
                        category = max(2, category)
        logger.debug('check exon %s resulted in %s', e, altsplice)
        return altsplice, category

    def _check_junction(self, j1, j2, e, e2):
        ''' check a junction in the segment graph

        * check presence e1-e2 junction in ref (-> if not exon skipping or novel junction)
            * presence is defined a direct junction from an ref exon (e.g. from self) overlapping e1 to an ref exon overlapping e2
        * AND find j3 and j4: first node overlapping e2  and last node overlapping e2
        * more specifically:
            * j3: first node ending after e2 start, or len(self)
            * j4: last node starting before e2 end (assuming there is such a node)'''
        altsplice = {}
        j3 = next((j for j in range(j2 + 1, len(self)) if self[j][1] > e2[0]), len(self))
        j4 = next((j - 1 for j in range(j3, len(self)) if self[j].start >= e2[1]), len(self) - 1)
        if j3 == len(self) or self[j3].start > e2[1]:
            return j3, j4, altsplice  # no overlap with e2
        if e[1] == self[j2].end and e2[0] == self[j3].start and j3 in self[j2].suc.values():
            return j3, j4, altsplice  # found direct junction
        # find skipped exons within e1-e2 intron
        exon_skipping = set()
        for j in range(j1, j2 + 1):
            for tr, j_suc in self[j].suc.items():
                if j_suc <= j2 or j_suc > j4:  # befor junction or exceding it
                    continue  # -> no exon skipping
                if self._pas[tr] < j4:  # transcripts ends within e1-e2 intron
                    continue  # -> no exon skipping
                j_spliced = self._get_next_spliced(tr, j)
                if j_spliced is None:  # end of transcript
                    continue  # -> no exon skipping
                if j_spliced < j3:  # found splice junction from e1 into e1-e2 intron
                    j_spliced_end = self._get_exon_end(tr, j_spliced)
                    if j_spliced_end < j3:  # ref exon ends within e1-e2 intron
                        exon_skipping.add(tr)  # found exon skipping
        if exon_skipping:  # path from e1 over another exon e_skip to e2 present
            # now we need to find the exon boundaries of the skipped exon
            exons = []
            e_start = e[1] - 1
            e_end = e[1]
            for j in range(j2 + 1, j3 + 1):
                if e_start and self[j].start > self[j - 1].end:  # intron befor j
                    exons.append([e_start, e_end])
                    e_start = 0
                if exon_skipping.intersection(self[j].suc):  # exon at j
                    if e_start == 0:  # new exon start
                        e_start = self[j].start
                    e_end = self[j].end
                elif e_start:  # intron at j
                    exons.append([e_start, e_end])
                    e_start = 0
            if len(exons) > 1:
                altsplice.setdefault('exon skipping', []).extend(exons[1:])
        elif e[1] == self[j2].end and e2[0] == self[j3].start:  # e1-e2 path is not present, but splice sites are
            altsplice.setdefault('novel junction', []).append([e[1], e2[0]])  # for example mutually exclusive exons spliced togeter

        logger.debug('check junction %s - %s resulted in %s', e[0], e[1], altsplice)

        return j3, j4, altsplice

    def fuzzy_junction(self, exons, size):
        '''Looks for "fuzzy junctions" in the provided transcript.

        For each intron from "exons", look for introns in the splice graph shifted by less than "size".
        These shifts may be produced by ambigious alignments.

        :param exons: A list of exon tuples representing the transcript
        :type exons: list
        :param size: The maximum size of the fuzzy junction
        :type size: int
        :return: a dict with the intron number as key and the shift as value (assuming size is smaller than introns)
        :rtype: dict'''
        fuzzy = {}
        if size < 1:  # no need to check
            return fuzzy
        j1 = 0
        idx = range(j1, len(self))
        for i, e1 in enumerate(exons[:-1]):
            e2 = exons[i + 1]
            try:
                # find j1: first node intersecting size range of e1 end
                if self[j1].end + min(size, e1[1] - e1[0]) < e1[1]:
                    j1 = next(j for j in idx if self[j].end + size >= e1[1])
            except (StopIteration, IndexError):  # transcript end - we are done
                break
            shift = []
            while j1 < len(self) and self[j1].end - e1[1] <= min(size, e1[1] - e1[0]):  # in case there are several nodes starting in the range around e1
                shift_e1 = self[j1].end - e1[1]
                # print(f'{i} {e1[1]}-{e2[0]} {shift_e1}')
                if shift_e1 == 0:  # no shift required at this intron
                    break
                if any(self[j2].start - e2[0] == shift_e1 for j2 in set(self[j1].suc.values())):
                    shift.append(shift_e1)
                j1 += 1
            else:  # junction not found in sg
                if shift:  # but shifted juction is present
                    fuzzy[i] = sorted(shift, key=abs)[0]  # if there are several possible shifts, provide the smallest
        return fuzzy

    def find_splice_sites(self, exons):
        '''Checks whether the splice sites of a new transcript are present in the segment graph.

        :param exons: A list of exon tuples representing the transcript
        :type exons: list
        :return: boolean array indicating whether the splice site is contained or not'''

        j = 0
        sites = np.zeros((len(exons) - 1) * 2, dtype=bool)
        for i, (e1, e2) in enumerate(pairwise(exons)):
            while(self[j].end < e1[1]):
                j += 1
                if j == len(self):
                    return sites
            if self[j].end == e1[1]:
                sites[i * 2] = True
            while(self[j].start < e2[0]):
                j += 1
                if j == len(self):
                    return sites
            if self[j].start == e2[0]:
                sites[i * 2 + 1] = True
        return sites

    def get_overlap(self, exons):
        '''Compute the exonic overlap of a new transcript with the segment graph.

        :param exons: A list of exon tuples representing the transcript
        :type exons: list
        :return: a tuple: the overlap with the gene, and a list of the overlaps with the transcripts'''

        ol = 0
        j = 0
        tr_ol = [0 for _ in self._pas]
        for e in exons:
            while(self[j].end < e[0]):  # no overlap, go on
                j += 1
                if j == len(self):
                    return ol, tr_ol
            while(self[j].start < e[1]):
                i_end = min(e[1], self[j].end)
                i_start = max(e[0], self[j].start)
                ol += (i_end - i_start)
                for trid in self[j].suc.keys():
                    tr_ol[trid] += (i_end - i_start)
                for trid, pas in enumerate(self._pas):
                    if pas == j:
                        tr_ol[trid] += (i_end - i_start)
                if self[j].end > e[1]:
                    break
                j += 1
                if j == len(self):
                    return ol, tr_ol

        return ol, tr_ol

    def get_intron_support_matrix(self, exons):
        '''Check the intron support for the provided transcript w.r.t. transcripts from self.

        This is supposed to be helpfull for the analysis of novel combinations of known splice sites.


        :param exons: A list of exon positions defining the transcript to check.
        :return: A boolean array of shape (n_transcripts in self)x(len(exons)-1).
            An entry is True iff the intron from "exons" is present in the respective transcript of self.'''
        node_iter = iter(self)
        ism = np.zeros((len(self._tss), len(exons) - 1), np.bool)  # the intron support matrix
        for intron_nr, (e1, e2) in enumerate(pairwise(exons)):
            try:
                node = next(n for n in node_iter if n.end >= e1[1])
            except StopIteration:
                return ism
            if node.end == e1[1]:
                for trid, suc in node.suc.items():
                    if self[suc].start == e2[0]:
                        ism[trid, intron_nr] = True
        return ism

    def get_exon_support_matrix(self, exons):
        '''Check the exon support for the provided transcript w.r.t. transcripts from self.

        This is supposed to be helpfull for the analysis of novel combinations of known splice sites.


        :param exons: A list of exon positions defining the transcript to check.
        :return: A boolean array of shape (n_transcripts in self)x(len(exons)-1).
            An entry is True iff the exon from "exons" is fully covered in the respective transcript of self.
            First and last exon are checked to overlap the first and last exon of the ref transcript but do not need to be fully covered'''
        esm = np.zeros((len(self._tss), len(exons)), np.bool)  # the intron support matrix
        for tr_nr, tss in enumerate(self._tss):  # check overlap of first exon
            for j in range(tss, len(self)):
                if overlap(self[j], exons[0]):
                    esm[tr_nr, 0] = True
                elif self[j].suc.get(tr_nr, None) == j + 1 and j - 1 < len(self) and self[j].end == self[j + 1].start:
                    continue
                break
        for tr_nr, pas in enumerate(self._pas):  # check overlap of last exon
            for j in range(pas, -1, -1):
                if overlap(self[j], exons[-1]):
                    esm[tr_nr, -1] = True
                elif self[j].pre.get(tr_nr, None) == j - 1 and j > 0 and self[j].start == self[j - 1].end:
                    continue
                break

        j2 = 0
        for e_nr, e in enumerate(exons[1:-1]):
            j1 = next((j for j in range(j2, len(self)) if self[j].end > e[0]))
            # j1: index of first segment ending after exon start (i.e. first overlapping segment)
            j2 = next((j - 1 for j in range(j1, len(self)) if self[j].start >= e[1]), len(self) - 1)
            # j2: index of last segment starting befor exon end (i.e. last overlapping segment)
            if self[j1].start <= e[0] and self[j2].end >= e[1]:
                covered = set.intersection(*(set(self[j].suc) for j in range(j1, j2 + 1)))
                if covered:
                    esm[covered, e_nr + 1] = True
        return esm

    def get_intersects(self, exons):
        '''Computes the splice junction exonic overlap of a new transcript with the segment graph.

        :param exons: A list of exon tuples representing the transcript
        :type exons: list
        :return: the splice junction overlap and exonic overlap'''

        intersect = [0, 0]
        i = j = 0
        while True:
            if self[j][0] == exons[i][0] and any(self[k][1] < self[j][0] for k in self[j].pre.values()):
                intersect[0] += 1  # same position and actual splice junction(not just tss or pas and internal junction)
            if self[j][1] == exons[i][1] and any(self[k][0] > self[j][1] for k in self[j].suc.values()):
                intersect[0] += 1
            if self[j][1] > exons[i][0] and exons[i][1] > self[j][0]:  # overlap
                intersect[1] += min(self[j][1], exons[i][1]) - max(self[j][0], exons[i][0])
            if exons[i][1] < self[j][1]:
                i += 1
            else:
                j += 1
            if i == len(exons) or j == len(self):
                return(intersect)

    @deprecated
    def _find_ts_candidates(self, coverage):
        '''Computes a metric indicating template switching.'''
        for i, gnode in enumerate(self._graph[:-1]):
            if self._graph[i + 1].start == gnode.end:  # jump candidates: introns that start within an exon
                jumps = {idx: n for idx, n in gnode.suc.items() if n > i + 1 and self._graph[n].start == self._graph[n - 1].end}
                # find jumps (n>i+1) and check wether they end within an exon begin(jumptarget)==end(node before)
                jump_weight = {}
                for idx, target in jumps.items():
                    jump_weight.setdefault(target, [0, []])
                    jump_weight[target][0] += coverage[:, idx].sum(0)
                    jump_weight[target][1].append(idx)

                for target, (w, idx) in jump_weight.items():
                    long_idx = set(idx for idx, n in gnode.suc.items() if n == i + 1) & set(idx for idx, n in self[target].pre.items() if n == target - 1)
                    try:
                        longer_weight = coverage[:, list(long_idx)].sum()
                    except IndexError:
                        print(long_idx)
                        raise
                    yield gnode.end, self[target].start, w, longer_weight, idx

    def _is_spliced(self, trid, ni1, ni2):
        'checks if transcript is spliced (e.g. has an intron) between nodes ni1 and ni2'
        if any(self[i].end < self[i + 1].start for i in range(ni1, ni2)):  # all transcripts are spliced
            return True
        if all(trid in self[i].suc for i in range(ni1, ni2)):
            return False
        return True

    def _get_next_spliced(self, trid, node):
        'find the next spliced node for given transcript'
        while node != self._pas[trid]:
            try:
                next_node = self[node].suc[trid]  # raises error if trid not in node
            except KeyError:
                logger.error('trid %s seems to be not in node %s', trid, node)
                raise
            if self[next_node].start > self[node].end:
                return next_node
            node = next_node
        return None

    def _get_exon_end(self, trid, node):
        'find the end of the exon to which node belongs for given transcript'
        while node != self._pas[trid]:
            try:
                next_node = self[node].suc[trid]  # raises error if trid not in node
            except KeyError:
                logger.error('trid %s seems to be not in node %s', trid, node)
                raise
            if self[next_node].start > self[node].end:
                return node
            node = next_node
        return node

    def _get_exon_start(self, trid, node):
        'find the start of the exon to which node belongs for given transcript'
        while node != self._tss[trid]:
            try:
                next_node = self[node].pre[trid]  # raises error if trid not in node
            except KeyError:
                logger.error('trid %s seems to be not in node %s', trid, node)
                raise
            if self[next_node].end < self[node].start:
                return node
            node = next_node
        return node

    def _find_splice_bubbles_at_position(self, tids, pos):
        '''function to refind bubbles at a certain genomic position.
        This turns out to be fundamentally different compared to iterating over all bubbles, hence it is a complete rewrite of the function.
        On the positive site, the functions can validate each other. I tried to reuse the variable names.
        If both functions yield same results, there is a good chance that the compex code is actually right.'''

        if any(t < 5 for t in tids):

            try:
                i, nA = next((idx, n) for idx, n in enumerate(self) if n.end == pos[0])
                if len(pos) == 3:
                    middle = [next(idx + i for idx, n in enumerate(self[i:]) if n.start > pos[1])]
                    j, nB = next((idx + middle[0], n) for idx, n in enumerate(self[middle[0]:]) if n.start == pos[2])
                else:
                    j, nB = next((idx + i, n) for idx, n in enumerate(self[i:]) if n.start == pos[1])
                    middle = [idx for idx in range(i + 2, j)]
            except StopIteration as e:
                raise ValueError(f"cannot find segments at {pos} in segment graph") from e

            direct = set()  # primary
            indirect = set(), set(), set(), set()  # for es, as at start, as at end, ir
            for tr, node_id in nA.suc.items():
                type_id = 0
                if tr not in nB.pre:
                    continue
                if node_id == j:
                    direct.add(tr)
                    continue
                if self[node_id].start == nA.end:
                    type_id += 2
                if self[nB.pre[tr]].end == nB.start:
                    type_id += 1
                indirect[type_id].add(tr)
            for t in tids:
                if t < 4 and direct and indirect[t]:  # ES,3AS.5AS,IR
                    yield direct, indirect[t], i, j, t
                elif t == 4 and len(indirect[0]) > 2:  # ME
                    me = list()
                    seen_alt = set()
                    for middle_idx in middle:
                        alt, prim = set(), set()
                        for tr in indirect[0]:  # spliced at nA and nB
                            if nB.pre[tr] < middle_idx:
                                alt.add(tr)
                            elif nA.suc[tr] >= middle_idx:
                                prim.add(tr)
                        if prim and alt - seen_alt:  # make sure there is at least one new alt transcript with this middle node.
                            me.append((prim, alt, i, j, t))
                            seen_alt.update(alt)
                    seen_prim = set()
                    for me_event in reversed(me):
                        if me_event[0] - seen_prim:  # report only if there is a new transcript in prim with respect to middle nodes to the right
                            yield me_event
                            seen_prim.update(me_event[0])
        if any(t > 4 for t in tids):
            try:
                i, nA = next((idx, n) for idx, n in enumerate(self) if n.start >= pos[0])
                j, nB = next(((idx + i, n) for idx, n in enumerate(self[i:]) if n.end >= pos[-1]), (len(self) - 1, self[-1]))
            except StopIteration as e:
                raise ValueError(f"cannot find segments at {pos} in segment graph") from e

            if 5 in tids and nB.end == pos[-1]:  # TSS on +, PAS on -
                alt = {tr for tr, tss in enumerate(self._tss) if i <= tss <= j and self._get_exon_end(tr, tss) == j}
                if alt:  # find compatible alternatives: end after tss /start before pas
                    prim = [tr for tr, pas in enumerate(self._pas) if tr not in alt and pas > j]  # prim={tr for tr in range(len(self._tss)) if tr not in alt}
                    if prim:
                        yield prim, alt, i, j, 5
            if 6 in tids and nA.start == pos[0]:
                alt = {tr for tr, pas in enumerate(self._pas) if i <= pas <= j and self._get_exon_start(tr, pas) == i}
                if alt:
                    prim = [tr for tr, tss in enumerate(self._tss) if tr not in alt and tss < i]
                    if prim:
                        yield prim, alt, i, j, 6

    def find_splice_bubbles(self, types=None, pos=None):
        '''Searches for alternative paths in the segment graph ("bubbles").

        Bubbles are defined as combinations of nodes x_s and x_e with more than one path from x_s to x_e.

        :param types: A tuple with event types to find. Valid types are ('ES','3AS', '5AS','IR' or 'ME', 'TSS', 'PAS').
            If ommited, all types are considered
        :param pos: If specified, restrict the search on specific position.
            This is useful to find the supporting transcripts for a given type if the position is known.

        :return: Tuple with 1) transcript indices of primary (e.g. most direct) paths and 2) alternative paths respectivly,
            as well as 3) start and 4) end node ids and 5) type of alternative event
            ('ES','3AS', '5AS','IR' or 'ME', 'TSS', 'PAS')'''

        if types is None:
            types = ('ES', '3AS', '5AS', 'IR', 'ME', 'TSS', 'PAS')
        elif isinstance(types, str):
            types = (types,)
        alt_types = ('ES', '5AS', '3AS', 'IR', 'ME', 'PAS', "TSS") if self.strand == '-' else ('ES', '3AS', '5AS', 'IR', 'ME', 'TSS', "PAS")

        if pos is not None:
            for prim, alt, i, j, alt_tid in self._find_splice_bubbles_at_position([i for i, t in enumerate(alt_types) if t in types], pos):
                yield list(prim), list(alt), i, j, alt_types[alt_tid]
            return
        if any(t in types for t in ('ES', '3AS', '5AS', 'IR', 'ME')):
            # alternative types: intron retention, alternative splice site at left and right, exon skipping, mutually exclusive
            inB_sets = [(set(), set())]  # list of spliced and unspliced transcripts joining in B
            # node_matrix=self.get_node_matrix()
            for i, nB in enumerate(self[1:]):
                inB_sets.append((set(), set()))
                unspliced = self[i].end == nB.start
                for tr, node_id in nB.pre.items():
                    inB_sets[i + 1][unspliced and node_id == i].add(tr)
            for i, nA in enumerate(self):
                junctions = sorted(list(set(nA.suc.values())))  # target nodes for junctions from node A ordered by intron size
                if len(junctions) < 2:
                    continue  # no alternative
                outA_sets = {}  # transcripts supporting the different  junctions
                for tr, node_id in nA.suc.items():
                    outA_sets.setdefault(node_id, set()).add(tr)
                unspliced = nA.end == self[junctions[0]].start
                alternative = ([], outA_sets[junctions[0]]) if unspliced else (outA_sets[junctions[0]], [])
                # nC_dict aims to avoid recalculation of nC for ME events
                # tr -> node at start of 2nd exon C for tr such that there is one exon (B) (and both flanking introns) between nA and C; None if transcript ends
                nC_dict = {}
                me_alt_seen = set()  # ensure that only ME events with novel tr are reported
                logger.debug('checking node %s: %s (%s)', i, nA, list(zip(junctions, [outA_sets[j] for j in junctions])))
                for j_idx, joi in enumerate(junctions[1:]):  # start from second, as first does not have an alternative
                    alternative = [{tr for tr in alternative[i] if self._pas[tr] > joi} for i in range(2)]  # check that transcripts extend beyond nB
                    logger.debug(alternative)
                    found = [trL1.intersection(trL2) for trL1 in alternative for trL2 in inB_sets[joi]]  # alternative transcript sets for the 4 types
                    #  5th type: mutually exclusive (outdated handling of ME for reference)
                    # found.append(set.union(*alternative)-inB_sets[joi][0]-inB_sets[joi][1])
                    logger.debug('checking junction %s (tr=%s) and found %s at B=%s', joi, outA_sets[joi], found, inB_sets[joi])
                    for alt_tid, alt in enumerate(found):
                        if alt_types[alt_tid] in types and alt:
                            yield list(outA_sets[joi]), list(alt), i, joi, alt_types[alt_tid]
                    # me_alt=set.union(*alternative)-inB_sets[joi][0]-inB_sets[joi][1] #search 5th type: mutually exclusive
                    if 'ME' in types:
                        me_alt = alternative[0] - inB_sets[joi][0] - inB_sets[joi][1]  # search 5th type: mutually exclusive - needs to be spliced
                        if me_alt - me_alt_seen:  # there is at least one novel alternative transcript
                            # for ME we need to find (potentially more than one) node C where the alternatives rejoin
                            for tr in me_alt:  # find node C for all me_alt
                                nC_dict.setdefault(tr, self._get_next_spliced(tr, nA.suc[tr]))
                                if nC_dict[tr] is None:  # transcript end in nB, no node C
                                    me_alt_seen.add(tr)  # those are not of interest for ME
                            inC_sets = {}  # dict of nC indices to sets of nA-intron-nB-intron-nC transcripts, where nB >= joi
                            for nB_i in junctions[j_idx + 1:]:  # all primary junctions
                                for tr in outA_sets[nB_i]:  # primary transcripts
                                    nC_dict.setdefault(tr, self._get_next_spliced(tr, nB_i))  # find nC
                                    if nC_dict[tr] is None:
                                        continue
                                    if nB_i == joi:  # first, all primary tr/nCs from joi added
                                        inC_sets.setdefault(nC_dict[tr], set()).add(tr)
                                    elif nC_dict[tr] in inC_sets:  # then add primary tr that also rejoin at any of the joi nC
                                        inC_sets[nC_dict[tr]].add(tr)
                                if not inC_sets:  # no nC for any of the joi primary tr - no need to check the other primaries
                                    break
                            for nC_i, me_prim in sorted(inC_sets.items()):
                                found_alt = {tr for tr in me_alt if nC_dict[tr] == nC_i}
                                if found_alt - me_alt_seen:  # ensure, there is a new alternative
                                    yield (list(me_prim), list(found_alt), i, nC_i, 'ME')
                                    # me_alt=me_alt-found_alt
                                    me_alt_seen.update(found_alt)
                    alternative[0].update(outA_sets[joi])  # now transcripts supporting joi join the alternatives
        if "TSS" in types or "PAS" in types:
            yield from self._find_start_end_events(types)

    def _find_start_end_events(self, types):
        '''Searches for alternative TSS/PAS in the segment graph.

        All transcripts sharing the same first/ last splice junction are considered to start/end at the same site and are summarized.
        Alternative transcripts are all other transcripts of the gene that end after the TSS or respectivly start before the PAS.

        :return: Tuple with 1) transcript ids sharing common start exon and 2) alternative transcript ids respectivly,
            as well as 3) start and 4) end node ids of the exon and 5) type of alternative event ("TSS" or "PAS")'''
        tss = {}
        pas = {}
        tss_start = {}
        pas_end = {}
        for tr, (start1, end1) in enumerate(zip(self._tss, self._pas)):
            start2 = self._get_exon_end(tr, start1)
            if start2 != end1:  # skip single exon transcripts
                # tss:  key: last node idx of start exon (e.g. first real splice junction),
                #       value: a list of transcripts sharing this node as the last node of their start exon)
                tss.setdefault(start2, set()).add(tr)
                # tss_start is first node of start exon (wrt all transcripts sharing same last node of start exon)
                tss_start[start2] = min(start1, tss_start.get(start2, start1))
                end2 = self._get_exon_start(tr, end1)
                pas.setdefault(end2, set()).add(tr)
                pas_end[end2] = max(end1, pas_end.get(end2, end1))
        alt_types = ['PAS', 'TSS'] if self.strand == '-' else ['TSS', 'PAS']
        if alt_types[0] in types:
            for n, tr_set in tss.items():  # find compatible alternatives: end after tss /start before pas
                alt_tr = [tr for tr, pas in enumerate(self._pas) if tr not in tr_set and pas > n]
                yield (alt_tr, list(tr_set), tss_start[n], n, alt_types[0])
        if alt_types[1] in types:
            for n, tr_set in pas.items():
                alt_tr = [tr for tr, tss in enumerate(self._tss) if tr not in tr_set and tss < n]
                yield (alt_tr, list(tr_set), n, pas_end[n], alt_types[1])

    def is_exonic(self, position):
        '''Checks whether the position is within an exon.

        :param position: The genomic position to check.
        :return: True, if the position overlaps with an exon, else False.'''
        for n in self:
            if n[0] <= position and n[1] >= position:
                return True
        return False

    def _get_all_exons(self, nodeX, nodeY, tr):
        'get all exons from nodeX to nodeY for transcripts tr'
        # TODO: add option to extend first and last exons
        n = max(nodeX, self._tss[tr])  # if tss>nodeX start there
        if tr not in self[n].pre and self._tss[tr] != n:  # nodeX is not part of tr
            for i in range(n, nodeY + 1):
                if tr in self[n].suc:
                    n = i
                    break
            else:
                return []
        if n > nodeY:
            return []
        exons = [[self[n].start, self[n].end]]
        while n < nodeY:
            try:
                n = self[n].suc[tr]
            except KeyError:  # end of transcript before nodeY was reached
                break
            if self[n].start == exons[-1][1]:
                exons[-1][1] = self[n].end
            else:
                exons.append([self[n].start, self[n].end])
        return [tuple(e) for e in exons]

    def __getitem__(self, key):
        return self._graph[key]

    def __len__(self):
        return len(self._graph)


class SegGraphNode(tuple):
    '''A node in a segment graph represents an exonic segment.'''
    def __new__(cls, start, end, pre=None, suc=None):
        if pre is None:
            pre = dict()
        if suc is None:
            suc = dict()
        return super(SegGraphNode, cls).__new__(cls, (start, end, pre, suc))

    def __getnewargs__(self):
        return tuple(self)

    @property
    def start(self) -> int:
        '''the (genomic 5') start of the segment'''
        return self.__getitem__(0)

    @property
    def end(self) -> int:
        '''the (genomic 3') end of the segment'''
        return self.__getitem__(1)

    @property
    def pre(self) -> dict:
        '''the predecessor segments of the segment (linked nodes downstream)'''
        return self.__getitem__(2)

    @property
    def suc(self) -> dict:
        '''the successor  segments of the segment (linked nodes downstream)'''
        return self.__getitem__(3)


class SpliceGraph():
    '''(Experimental) Splice Graph Implementation

    Nodes represent splice sites and are tuples of genomic positions and a "lower" flag.
    The "lower flag" is true, if the splice site is a genomic 5' end of an exon
    Nodes are kept sorted, so iteration over splicegraph returns all nodes in genomic order
    Edges are assessed with SpliceGraph.pre(node, [tr_nr]) and SpliceGraph.suc(node, [tr_nr]) functions.
    If no tr_nr is provided, a dict with all incoming/outgoing edges is returned'''
    # @experimental

    def __init__(self, is_reverse, graph, fwd_starts, rev_starts):
        self.is_reverse = is_reverse
        self._graph = graph
        self._fwd_starts = fwd_starts
        self._rev_starts = rev_starts

    @classmethod
    def from_transcript_list(cls, transcripts, strand):
        '''Compute the splice graph from a list of transcripts

        :param transcripts: A list of transcripts, which are lists of exons, which in turn are (start,end) tuples
        :type transcripts: list
        :param strand: the strand of the gene, either "+" or "-"
        :return: The SpliceGraph object
        :rtype: SpliceGraph'''

        assert strand in '+-', 'strand must be either "+" or "-"'

        graph = SortedDict()
        fwd_starts = [(tr[0][0], True) for tr in transcripts]  # genomic 5'
        rev_starts = [(tr[-1][1], False) for tr in transcripts]  # genomic 3'

        for tr_nr, tr in enumerate(transcripts):
            graph.setdefault((tr[0][0], True), ({}, {}))

            for i, (b1, b2) in enumerate(pairwise(pos for e in tr for pos in e)):
                graph.setdefault((b2, bool(i % 2)), ({}, {}))
                graph[b2, bool(i % 2)][1][tr_nr] = b1, not bool((i) % 2)  # successor
                graph[b1, not bool(i % 2)][0][tr_nr] = b2, bool((1) % 2)  # predesessor
        sg = cls(strand == '-', graph, fwd_starts, rev_starts)
        return sg

    def __iter__(self):
        return self._graph.__iter__()

    def __len__(self):
        return self._graph.__len__()

    @experimental
    def add(self, tr) -> None:
        '''Add one transcript to the existing graph.

        :param tr: A list of exon tuples representing the transcript to add
        :type tr: list'''

        tr_nr = len(self._fwd_starts)
        self._fwd_starts.append((tr[0][0], True))  # genomic 5'
        self._rev_starts.append((tr[-1][1], False))  # genomic 3'
        self._graph.setdefault((tr[0][0], True), ({}, {}))
        for i, (b1, b2) in enumerate(pairwise(pos for e in tr for pos in e)):
            self._graph.setdefault((b2, bool(i % 2)), ({}, {}))
            self._graph[b2, bool(i % 2)][1][tr_nr] = b1, not bool((i) % 2)  # successor
            self._graph[b1, not bool(i % 2)][0][tr_nr] = b2, bool((1) % 2)

    def suc(self, node, tr_nr=None) -> Union[int, dict]:
        '''get index of successor node (next genomic upstream node) of transcript, or, if tr_nr is omitted, a dict with successors for all transcripts.

        :param node: index of the originating node
        :type node: int
        :param tr_nr: index of the transcript (optional)
        :type tr_nr: int'''
        edges = self._graph[node][0]
        if tr_nr is None:
            return edges
        return edges[tr_nr]

    def pre(self, node, tr_nr=None) -> Union[int, dict]:
        '''get index of predesessor node (next genomic downstream node) of transcript, or, if tr_nr is omitted, a dict with predesessors for all transcripts.

        :param node: index of the originating node
        :type node: int
        :param tr_nr: index of the transcript (optional)
        :type tr_nr: int'''

        edges = self._graph[node][1]
        if tr_nr is None:
            return edges
        return edges[tr_nr]
