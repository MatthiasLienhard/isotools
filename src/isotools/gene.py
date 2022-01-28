from intervaltree import Interval
from collections.abc import Iterable
from Bio.Seq import Seq
import numpy as np
import copy
import itertools
from .splice_graph import SegmentGraph
from .short_read import Coverage
from ._transcriptome_filter import SPLICE_CATEGORY
from ._utils import pairwise, _filter_event
from ._transcriptome_stats import pairwise_event_test

import logging
logger = logging.getLogger('isotools')


class Gene(Interval):
    'This class stores all gene information and transcripts. It is derived from intervaltree.Interval.'
    required_infos = ['ID', 'chr', 'strand']

    # initialization
    def __new__(cls, begin, end, data, transcriptome):
        return super().__new__(cls, begin, end, data)  # required as Interval (and Gene) is immutable

    def __init__(self, begin, end, data, transcriptome):
        self._transcriptome = transcriptome

    def __str__(self):
        return('Gene {} {}({}), {} reference transcripts, {} expressed transcripts'.format(
            self.name, self.region, self.strand, self.n_ref_transcripts, self.n_transcripts))

    def __repr__(self):
        return object.__repr__(self)

    from ._gene_plots import sashimi_plot, gene_track, sashimi_plot_short_reads, sashimi_figure

    def short_reads(self, idx):
        '''Returns the short read coverage profile for a short read sample.

        :param idx: The index of the short read sample. '''

        try:
            return self.data['short_reads'][idx]
        except (KeyError, IndexError):
            srdf = self._transcriptome.infos['short_reads']  # raises key_error if no short reads added
            self.data.setdefault('short_reads', [])
            for i in range(len(self.data['short_reads']), len(srdf)):
                self.data['short_reads'].append(Coverage.from_bam(srdf.file[i], self))
        return self.data['short_reads'][idx]

    def correct_fuzzy_junctions(self, trid, size, modify=True):
        '''Corrects for splicing shifts.

         This function looks for "shifted junctions", e.g. same difference compared to reference annotaiton at both donor and acceptor)
         presumably caused by ambigous alignments. In these cases the positions are adapted to the reference position (if modify is set).

         :param trid: The index of the transcript to be checked.
         :param size: The maximum shift to be corrected.
         :param modify: If set, the exon positions are corrected according to the reference.'''

        exons = trid['exons']
        shifts = self.ref_segment_graph.fuzzy_junction(exons, size)
        if shifts and modify:
            for i, sh in shifts.items():
                if exons[i][0] <= exons[i][1] + sh and exons[i + 1][0] + sh <= exons[i + 1][1]:
                    exons[i][1] += sh
                    exons[i + 1][0] += sh
            trid['exons'] = [e for e in exons if e[0] < e[1]]  # remove zero length exons
        return shifts

    def _to_gtf(self, trids, source='isoseq'):
        '''Creates the gtf lines of the gene as strings.'''
        donotshow = {'transcripts', 'short_exons', 'segment_graph'}
        info = {'gene_id': self.id, 'gene_name': self.name}
        lines = [None]
        starts = []
        ends = []
        for i in trids:
            tr = self.transcripts[i]
            info['transcript_id'] = f'{info["gene_id"]}_{i}'
            starts.append(tr['exons'][0][0] + 1)
            ends.append(tr['exons'][-1][1])
            trinfo = info.copy()
            if 'downstream_A_content' in tr:
                trinfo['downstream_A_content'] = f'{tr["downstream_A_content"]:0.3f}'
            if tr['annotation'][0] == 0:  # FSM
                refinfo = {}
                for refid in tr['annotation'][1]['FSM']:
                    for k in self.ref_transcripts[refid]:
                        if k == 'exons':
                            continue
                        elif k == 'CDS':
                            if self.strand == '+':
                                cds_start, cds_end = self.ref_transcripts[refid]['CDS']
                            else:
                                cds_end, cds_start = self.ref_transcripts[refid]['CDS']
                            refinfo.setdefault('CDS_start', []).append(str(cds_start))
                            refinfo.setdefault('CDS_end', []).append(str(cds_end))
                        else:
                            refinfo.setdefault(k, []).append(str(self.ref_transcripts[refid][k]))
                for k, vlist in refinfo.items():
                    trinfo[f'ref_{k}'] = ','.join(vlist)
            else:
                trinfo['novelty'] = ','.join(k for k in tr['annotation'][1])
            lines.append((self.chrom, source, 'transcript', tr['exons'][0][0] + 1, tr['exons'][-1][1], '.',
                         self.strand, '.', '; '.join(f'{k} "{v}"' for k, v in trinfo.items())))
            noncanonical = tr.get('noncanonical_splicing', [])
            for enr, pos in enumerate(tr['exons']):
                exon_info = info.copy()
                exon_info['exon_id'] = f'{info["gene_id"]}_{i}_{enr}'
                if enr in noncanonical:
                    exon_info['noncanonical_donor'] = noncanonical[enr][:2]
                if enr+1 in noncanonical:
                    exon_info['noncanonical_acceptor'] = noncanonical[enr+1][2:]
                lines.append((self.chrom, source, 'exon', pos[0] + 1, pos[1], '.', self.strand, '.', '; '.join(f'{k} "{v}"' for k, v in exon_info.items())))
        if len(lines) > 1:
            # add gene line
            if 'reference' in self.data:
                info.update({k: v for k, v in self.data['reference'].items() if k not in donotshow})  # add reference gene specific fields
            lines[0] = (self.chrom, source, 'gene', min(starts), max(ends), '.', self.strand, '.', '; '.join(f'{k} "{v}"' for k, v in info.items()))
            return lines
        return []

    def add_noncanonical_splicing(self, genome_fh):
        '''Add information on noncanonical splicing.

        For all transcripts of the gene, scan for noncanonical (i.e. not GT-AG) splice sites.
        If noncanonical splice sites are present, the corresponding intron index (in genomic orientation) and the sequence
        i.e. the dinucleotides of donor and aceptor as XX-YY string are stored in the "noncannoncical_splicing" field of the transcript dicts.
        True noncanonical splicing is rare, thus it might indicate technical artifacts (template switching, missalignment, ...)

        :param genome_fh: A file handle of the genome fasta file.'''
        ss_seq = {}
        for tr in self.transcripts:
            pos = [(tr['exons'][i][1], tr['exons'][i + 1][0] - 2) for i in range(len(tr['exons']) - 1)]
            new_ss_seq = {site: genome_fh.fetch(self.chrom, site, site + 2).upper() for intron in pos for site in intron if site not in ss_seq}
            if new_ss_seq:
                ss_seq.update(new_ss_seq)

            if self.strand == '+':
                sj_seq = [ss_seq[d] + ss_seq[a] for d, a in pos]
            else:
                sj_seq = [str(Seq(ss_seq[d] + ss_seq[a]).reverse_complement()) for d, a in pos]

            nc = [(i, seq) for i, seq in enumerate(sj_seq) if seq != 'GTAG']
            if nc:
                tr['noncanonical_splicing'] = nc

    def add_direct_repeat_len(self, genome_fh, delta=15, max_mm=2, wobble=2):
        '''Computes direct repeat length.

        This function counts the number of consequtive equal bases at donor and acceptor sites of the splice junctions.
        This information is stored in the "direct_repeat_len" filed of the transcript dictionaries.
        Direct repeats longer than expected by chance indicate template switching.

        :param genome_fh: The file handle to the genome fasta.
        :param delta: The maximum length of direct repeats that can be found.
        :param max_mm: The maximum length of direct repeats that can be found.
        :param wobble: The maximum length of direct repeats that can be found.'''

        intron_seq = {}
        score = {}

        for tr in self.transcripts:
            for intron in ((tr['exons'][i][1], tr['exons'][i + 1][0]) for i in range(len(tr['exons']) - 1)):
                for pos in intron:
                    try:
                        intron_seq.setdefault(pos, genome_fh.fetch(self.chrom, pos - delta, pos + delta))
                    except (ValueError, IndexError):  # N padding at start/end of the chromosomes
                        chr_len = genome_fh.get_reference_length(self.chrom)
                        seq = genome_fh.fetch(self.chrom, max(0, pos - delta), min(chr_len, pos + delta))
                        if pos - delta < 0:
                            seq = ''.join(['N'] * (pos - delta)) + seq
                        if pos + delta > chr_len:
                            seq += ''.join(['N'] * (pos + delta - chr_len))
                        intron_seq.setdefault(pos, seq)
                if intron not in score:
                    score[intron] = repeat_len(intron_seq[intron[0]], intron_seq[intron[1]], wobble=wobble, max_mm=max_mm)

        for tr in self.transcripts:
            tr['direct_repeat_len'] = [min(score[(e1[1], e2[0])], delta) for e1, e2 in pairwise(tr['exons'])]

    def add_threeprime_a_content(self, genome_fh, length=30):
        '''Adds the information of the genomic A content downstream the transcript.

        High values of genomic A content indicate internal priming and hence genomic origin of the LRTS read.
        This function populates the 'downstream_A_content' field of the transcript dictionaries.

        :param geneome_fh: A file handle for the indexed genome fasta file.
        :param length: The length of the downstream region to be considered.
        '''
        a_content = {}
        for tr in (t for tL in (self.transcripts, self.ref_transcripts) for t in tL):
            if self.strand == '+':
                pos = tr['exons'][-1][1]
            else:
                pos = tr['exons'][0][0] - length
            if pos not in a_content:
                seq = genome_fh.fetch(self.chrom, max(0, pos), pos + length)
                if self.strand == '+':
                    a_content[pos] = seq.upper().count('A') / length
                else:
                    a_content[pos] = seq.upper().count('T') / length
            tr['downstream_A_content'] = a_content[pos]

    def add_fragments(self):
        '''Checks for transcripts that are fully contained in other transcripts.

        Transcripts that are fully contained in other transcripts are potential truncations.
        This function populates the 'fragment' filed of the transcript dictionaries with the indices of the containing transcripts,
        and the exon ids that match the first and last exons.'''

        for trid, containers in self.segment_graph.find_fragments().items():
            self.transcripts[trid]['fragments'] = containers  # list of (containing transcript id, first 5' exons, first 3'exons)

    def coding_len(self, trid):
        '''Returns length of 5\'UTR, coding sequence and 3\'UTR.

        :param trid: The transcript index for which the coding length is requested. '''

        try:
            exons = self.transcripts[trid]['exons']
            cds = self.transcripts[trid]['CDS']
        except KeyError:
            return None
        else:
            coding_len = _coding_len(exons, cds)
        if self.strand == '-':
            coding_len.reverse()
        return coding_len

    def get_infos(self, trid, keys, sample_i, group_i, **kwargs):
        '''Returns the transcript information specified in "keys" as a list.'''
        return [value for k in keys for value in self._get_info(trid, k, sample_i, group_i)]

    def _get_info(self, trid, key, sample_i, group_i, **kwargs):
        # returns tuples (as some keys return multiple values)
        if key == 'length':
            return sum((e - b for b, e in self.transcripts[trid]['exons'])),
        elif key == 'n_exons':
            return len(self.transcripts[trid]['exons']),
        elif key == 'exon_starts':
            return ','.join(str(e[0]) for e in self.transcripts[trid]['exons']),
        elif key == 'exon_ends':
            return ','.join(str(e[1]) for e in self.transcripts[trid]['exons']),
        elif key == 'annotation':
            # sel=['sj_i','base_i', 'as']
            if 'annotation' not in self.transcripts[trid]:
                return ('NA',) * 2
            nov_class, subcat = self.transcripts[trid]['annotation']
            # subcat_string = ';'.join(k if v is None else '{}:{}'.format(k, v) for k, v in subcat.items())
            return SPLICE_CATEGORY[nov_class], ','.join(subcat)  # only the names of the subcategories
        elif key == 'coverage':
            return self.coverage[sample_i, trid]
        elif key == 'tpm':
            return self.tpm(kwargs.get('pseudocount', 1))[sample_i, trid]
        elif key == 'group_coverage_sum':
            return tuple(self.coverage[si, trid].sum() for si in group_i)
        elif key == 'group_tpm_mean':
            return tuple(self.tpm(kwargs.get('pseudocount', 1))[si, trid].mean() for si in group_i)
        elif key in self.transcripts[trid]:
            val = self.transcripts[trid][key]
            if isinstance(val, Iterable):  # iterables get converted to string
                return str(val),
            else:
                return val,  # atomic (e.g. numeric)
        return 'NA',

    def _set_coverage(self, force=False):
        samples = self._transcriptome.samples
        cov = np.zeros((len(samples), self.n_transcripts), dtype=int)
        if not force:  # keep the segment graph if no new transcripts
            known = self.data.get('coverage', None)
            if known is not None and known.shape[1] == self.n_transcripts:
                if known.shape == cov.shape:
                    return
                cov[:known.shape[0], :] = known
                for i in range(known.shape[0], len(samples)):
                    for j, tr in enumerate(self.transcripts):
                        cov[i, j] = tr['coverage'].get(samples[i], 0)
                self.data['coverage'] = cov
                return
        for i, sa in enumerate(samples):
            for j, tr in enumerate(self.transcripts):
                cov[i, j] = tr['coverage'].get(sa, 0)
        self.data['coverage'] = cov
        self.data['segment_graph'] = None

    def tpm(self, pseudocount=1):
        '''Returns the transcripts per million (TPM).

        TPM is returned as a numpy array, with samples in columns and transcript isoforms in the rows.'''
        return (self.coverage+pseudocount)/self._transcriptome.sample_table['nonchimeric_reads'].values.reshape(-1, 1)*1e6

    @property
    def coverage(self):
        '''Returns the transcript coverage.

        Coverage is returned as a numpy array, with samples in columns and transcript isoforms in the rows.'''
        cov = self.data.get('coverage', None)
        if cov is not None:
            return cov
        self._set_coverage()
        return self.data['coverage']

    @property
    def gene_coverage(self):
        '''Returns the gene coverage.

        Total Coverage of the gene for each sample.'''

        return self.coverage.sum(1)

    @property
    def chrom(self):
        '''Returns the genes chromosome.'''
        return self.data['chr']

    @property
    def start(self):  # alias for begin
        return self.begin

    @property
    def region(self):
        '''Returns the region of the gene as a string in the format "chr:start-end".'''
        try:
            return '{}:{}-{}'.format(self.chrom, self.start, self.end)
        except KeyError:
            raise

    @property
    def id(self):
        '''Returns the gene id'''
        try:
            return self.data['ID']
        except KeyError:
            logger.error(self.data)
            raise

    @property
    def name(self):
        '''Returns the gene name'''
        try:
            return self.data['name']
        except KeyError:
            return self.id  # e.g. novel genes do not have a name (but id)

    @property
    def is_annotated(self):
        '''Returns "True" iff reference annotation is present for the gene.'''
        return 'reference' in self.data

    @property
    def is_expressed(self):
        '''Returns "True" iff gene is covered by at least one long read in at least one sample.'''
        return bool(self.transcripts)

    @property
    def strand(self):
        '''Returns the strand of the gene, e.g. "+" or "-"'''
        return self.data['strand']

    @property
    def transcripts(self):
        '''Returns the list of transcripts of the gene, as found by LRTS.'''
        try:
            return self.data['transcripts']
        except KeyError:
            return []

    @property
    def ref_transcripts(self):
        '''Returns the list of reference transcripts of the gene.'''
        try:
            return self.data['reference']['transcripts']
        except KeyError:
            return []

    @property
    def n_transcripts(self):
        '''Returns number of transcripts of the gene, as found by LRTS.'''
        return len(self.transcripts)

    @property
    def n_ref_transcripts(self):
        '''Returns number of reference transcripts of the gene.'''
        return len(self.ref_transcripts)

    @property
    def ref_segment_graph(self):  # raises key error if not self.is_annotated
        '''Returns the segment graph of the reference transcripts for the gene'''

        assert self.is_annotated, "reference segment graph requested on novel gene"
        if 'segment_graph' not in self.data['reference'] or self.data['reference']['segment_graph'] is None:
            exons = [tr['exons'] for tr in self.ref_transcripts]
            self.data['reference']['segment_graph'] = SegmentGraph(exons, self.strand)
        return self.data['reference']['segment_graph']

    @property
    def segment_graph(self):
        '''Returns the segment graph of the LRTS transcripts for the gene'''
        if 'segment_graph' not in self.data or self.data['segment_graph'] is None:
            exons = [tr['exons'] for tr in self.transcripts]
            self.data['segment_graph'] = SegmentGraph(exons, self.strand)
        return self.data['segment_graph']

    def __copy__(self):
        return Gene(self.start, self.end, self.data, self._transcriptome)

    def __deepcopy__(self, memo):  # does not copy _transcriptome!
        return Gene(self.start, self.end, copy.deepcopy(self.data, memo), self._transcriptome)

    def __reduce__(self):
        return Gene, (self.start, self.end, self.data, self._transcriptome)

    def copy(self):
        'Returns a shallow copy of self.'
        return self.__copy__()

    def gene_coordination_test(self, test="chi2", min_dist=1, min_total=100,
                               min_alt_fraction=.1, event_type=("ES", "5AS", "3AS", "IR", "ME")):
        '''
        Performs pairwise independence test for all pairs of Alternative Splicing Events (ASEs) in a gene.

        For all pairs of ASEs in a gene creates a contingency table and performs an indeppendence test.
        All ASEs A have two states, pri_A and alt_A, the primary and the alternative state respectivley.
        Thus, given two events A and B, we have four possible ways in which these events can occur on
        a transcript, that is, pri_A and pri_B, pri_A and alt_B, alt_A and pri_B, and alt_A and alt_B.
        These four values can be put in a contingency table and independence, or coordination,
        between the two events can be tested.

        :param test: Test to be performed. One of ("chi2", "fisher")
        :type test: str
        :param min_dist: Minimum distance (in nucleotides) between the two 
        Alternative Splicing Events for the pair to be tested
        :type min_dist: int
        :param min_total: The minimum total number of reads for an event to pass the filter
        :type min_total: int
        :param min_alt_fraction: The minimum fraction of read supporting the alternative
        :type min_alt_frction: float
        :param event_type:  A tuple with event types to test. Valid types are
        (‘ES’,’3AS’, ‘5AS’,’IR’ or ‘ME’, ‘TSS’, ‘PAS’). Default is ("ES", "5AS", "3AS", "IR", "ME")

        :return: a tuple (p_value, stat, Gene, ASE1_type, ASE2_type, ASE1_start, ASE1_end, ASE2_start,
        ASE2_end, priA_priB, priA_altB, altA_priB, altA_altB), where each element is a list
        containing the p_values, the statistics, the gene name, the type of the first ASE, 
        the type of the second ASE, the starting coordinate of the first ASE, the ending coordinate
        of the first ASE, the starting coordinate of the second ASE, the ending coordinate of the
        second ASE, and the four entries of the contingency table. 
        '''

        sg = self.segment_graph

        events = list(sg.find_splice_bubbles(types=event_type))
        events = [e for e in events if _filter_event(self, e, min_total=min_total,
                                                     min_alt_fraction=min_alt_fraction)]

        p_value = []
        stat = []
        Gene = []
        ASE1_type, ASE2_type = [], []
        ASE1_start = []
        ASE1_end = []
        ASE2_start = []
        ASE2_end = []
        priA_priB = []
        priA_altB = []
        altA_priB = []
        altA_altB = []

        for i, j in list(itertools.combinations(range(len(events)), 2)):

            test_res = pairwise_event_test(events[i], events[j], self, test=test, min_dist=min_dist, sg=sg)

            if test_res is not None:  # it is none for pairs of events whose distanca is smaller than min_dist
                if test_res[0] != 1:
                    p_value.append(test_res[0])
                    stat.append(test_res[1])
                    Gene.append(self.id)
                    ASE1_type.append(events[i][4])
                    ASE2_type.append(events[j][4])
                    ASE1_start.append(sg[events[i][2]].start)  # starting coordinate of event 1
                    ASE1_end.append(sg[events[i][3]].end)  # ending coordinate of event 1
                    ASE2_start.append(sg[events[j][2]].start)  # starting coordinate of event 2
                    ASE2_end.append(sg[events[j][3]].end)  # ending coordinate of event 2
                    priA_priB.append(test_res[2])
                    priA_altB.append(test_res[3])
                    altA_priB.append(test_res[4])
                    altA_altB.append(test_res[5])

        if len(p_value) == 0:
            return None

        return p_value, stat, Gene, ASE1_type, ASE2_type, ASE1_start, ASE1_end, ASE2_start, ASE2_end, priA_priB, priA_altB, altA_priB, altA_altB


def _coding_len(exons, cds):
    coding_len = [0, 0, 0]
    state = 0
    for e in exons:
        if state < 2 and e[1] >= cds[state]:
            coding_len[state] += cds[state] - e[0]
            if state == 0 and cds[1] <= e[1]:  # special case: CDS start and end in same exon
                coding_len[1] = cds[1] - cds[0]
                coding_len[2] = e[1] - cds[1]
                state += 2
            else:
                coding_len[state + 1] = e[1] - cds[state]
                state += 1
        else:
            coding_len[state] += e[1] - e[0]
    return coding_len


def repeat_len(seq1, seq2, wobble, max_mm):
    ''' Calcluate direct repeat length between seq1 and seq2
    '''
    score = [0]*(2*wobble+1)
    delta = int(len(seq1)/2-wobble)
    for w in range(2*wobble+1):  # wobble
        s1 = seq1[w:len(seq1)-(2*wobble-w)]
        s2 = seq2[wobble:len(seq2)-wobble]
        align = [a == b for a, b in zip(s1, s2)]
        score_left = find_runlength(reversed(align[:delta]), max_mm)
        score_right = find_runlength(align[delta:], max_mm)
        score[w] = max([score_left[fmm]+score_right[max_mm-fmm] for fmm in range(max_mm+1)])
    return max(score)


def find_runlength(align, max_mm):
    '''Find the runlength, e.g. the number of True in the list before the max_mm+1 False occur.
    '''
    score = [0]*(max_mm+1)
    mm = 0
    for a in align:
        if not a:
            mm += 1
            if mm > max_mm:
                return score
            score[mm] = score[mm-1]
        else:
            score[mm] += 1
    for i in range(mm+1, max_mm+1):
        score[i] = score[i-1]
    return score
