# methods and classes for the integration of short read data
from pysam import AlignmentFile
from ._utils import junctions_from_cigar
import numpy as np
import logging
logger = logging.getLogger('isotools')


class Coverage:
    'stores the illumina read coverage of a gene'
    # plan: make a binned version, or use run length encoding

    def __init__(self, cov, junctions, offset, chrom=None):
        self._cov = cov
        self._junctions = junctions
        self.reg = None if cov is None else (chrom, offset, offset+len(cov))
        self.bam_fn = None

    @classmethod
    def from_bam(cls, bam_fn, g, load=False):
        'assign the bam file'
        if load:
            with AlignmentFile(bam_fn, 'rb') as align:
                return cls.from_alignment(g, align)
        else:  # load on demand
            obj = cls.__new__(cls)
            obj._cov = None
            obj._junctions = None
            obj.bam_fn = bam_fn
            start = min(g.start, *[tr['exons'][0][0] for tr in g.transcripts])
            end = max(g.end, *[tr['exons'][-1][1] for tr in g.transcripts])
            obj.reg = (g.chrom, start, end)
            return obj

    @classmethod
    def from_alignment(cls, align_fh, g):
        'load the coverage from bam file'
        start = min(g.start, *[tr['exons'][0][0] for tr in g.transcripts])
        end = max(g.end, *[tr['exons'][-1][1] for tr in g.transcripts])
        cov, junctions = cls._import_coverage(align_fh, (g.chrom, start, end))
        obj = cls.__new__(cls)
        obj.__init__(cov, junctions, start)
        return obj

    @classmethod  # this is slow - called only if coverage is requested
    def _import_coverage(cls, align_fh, reg):
        delta = np.zeros(reg[2]-reg[1])
        junctions = {}
        for read in align_fh.fetch(*reg):
            exons = junctions_from_cigar(read.cigartuples, read.reference_start)
            # alternative: read.get_blocks() should be more efficient.. -todo: is it different?
            for i, exon in enumerate(exons):
                s = max(reg[1], min(reg[2]-1, exon[0]))-reg[1]
                e = max(reg[1], min(reg[2]-1, exon[1]))-reg[1]
                delta[s] += 1
                delta[e] -= 1
                if i > 0:
                    jpos = (exons[i-1][1], exon[0])
                    if jpos[1]-jpos[0] < 1 or jpos[0] < reg[1] or jpos[1] > reg[2]:
                        continue

                    junctions[jpos] = junctions.get(jpos, 0)+1
        cov = np.cumsum(delta)  # todo: use rle instead?
        return cov, junctions

    def load(self):
        'load the coverage from bam file'
        with AlignmentFile(self.bam_fn, 'rb') as align:
            logger.debug(f'Illumina coverage of region {self.reg[0]}:{self.reg[1]}-{self.reg[2]} is loaded from {self.bam_fn}')  # info or debug?
            self._cov, self._junctions = type(self)._import_coverage(align, self.reg)

    @property
    def junctions(self):
        if self._junctions is None:
            self.load()
        return self._junctions

    @property
    def profile(self):
        if self._cov is None:
            self.load()
        return self._cov  # todo: implement rle?

    def __getitem__(self, subscript):
        if isinstance(subscript, slice):
            return self.profile[slice(None if subscript.start is None else subscript.start-self.reg[1],
                                      None if subscript.stop is None else subscript.stop-self.reg[1],
                                      subscript.step)]  # does not get extended if outside range
        elif subscript < self.reg[1] or subscript >= self.reg[2]:
            logger.warning('requested coverage outside range')
            return None
        else:
            return(self.profile[subscript-self.reg[1]])
