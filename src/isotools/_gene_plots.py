
import collections.abc
import matplotlib.colors as plt_col
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
from math import log10
from ._utils import has_overlap, pairwise
import logging
logger = logging.getLogger('isotools')


def _label_overlap(pos1, pos2, width, height):
    if abs(pos1[0] - pos2[0]) < width and abs(pos1[1] - pos2[1]) < height:
        return True
    return False


DEFAULT_JPARAMS = [{'color': 'lightgrey', 'lwd': 1, 'draw_label': False},  # low coverage junctions
                   {'color': 'green', 'lwd': 1, 'draw_label': True},  # high coverage junctions
                   {'color': 'purple', 'lwd': 2, 'draw_label': True}]  # junctions of interest
DEFAULT_PARAMS = dict(min_cov_th=.001, high_cov_th=.05, text_width=.02, arc_type='both', text_height=1, exon_color='green')
DOMAIN_COLS = {"Family": "red", "Domain": "green", "Repeat": "orange", "Coiled-coil": "blue", "Motif": "grey", "Disordered": "pink"}


def extend_params(params):
    if params is None:
        params = dict()
    params.setdefault('jparams', [{}, {}, {}])
    # jparams=[params.pop(k,jparams[i]) for i,k in enumerate(['low_cov_junctions','high_cov_junctions','interest_junctions'])]
    for i, k1 in enumerate(['low_cov_junctions', 'high_cov_junctions', 'interest_junctions']):
        params['jparams'][i] = params.pop(k1, params['jparams'][i])
        for k2, v in DEFAULT_JPARAMS[i].items():
            params['jparams'][i].setdefault(k2, v)
    for k, v in DEFAULT_PARAMS.items():
        params.setdefault(k, v)
    return params


def get_index(samples, names):
    if not samples:
        return []
    if isinstance(names, list):
        idx = {sa: i for i, sa in enumerate(names)}
    else:
        idx = {sa: i for i, sa in names.items()}
    try:
        sidx = [idx[sa] for sa in samples]
    except KeyError:
        notfound = [sa for sa in samples if sa not in idx]
        logger.error('did not find the following samples: %s', ','.join(notfound))
        raise
    return sidx

# sashimi plots


def sashimi_figure(self, samples=None, short_read_samples=None, draw_gene_track=True, draw_other_genes=True,
                   long_read_params=None, short_read_params=None, junctions_of_interest=None, x_range=None):
    '''Arranges multiple Sashimi plots of the gene.

    The Sashimi figure consist of a reference gene track, long read sashimi plots for one or more samples or groups of samples,
    and optionally short read sashimi plots for one or more samples or groups of samples.

    :param samples: Definition of samples (as a list) or groups of samples (as a dict) for long read plots.
    :param short_read_samples: Definition of samples (as a list) or groups of samples (as a dict) for short read plots.
    :param draw_gene_track: Specify whether to plot the reference gene track.
    :param draw_other_genes: Specify, whether to draw other genes in the gene track.
    :param long_read_params: Dict with parameters for the long read plots, get passed to self.sashimi_plot.
        See isotools._gene_plots.DEFAULT_PARAMS and isotools._gene_plots.DEFAULT_JPARAMS
    :param short_read_params: Dict with parameters for the short read plots, get passed to self.sashimi_plot_short_reads.
        See isotools._gene_plots.DEFAULT_PARAMS and isotools._gene_plots.DEFAULT_JPARAMS
    :param junctions_of_interest: List of int pairs to define junctions of interest (which are highlighed in the plots)
    :param x_range: Genomic positions to specify the x range of the plot.
    :return: Tuple with figure and axses'''

    draw_gene_track = bool(draw_gene_track)

    if samples is None:
        samples = {}
    if short_read_samples is None:
        short_read_samples = {}
    if not samples and not short_read_samples:
        samples = {'all': None}
    if long_read_params is None:
        long_read_params = {}
    if short_read_params is None:
        short_read_params = {}

    f, axes = plt.subplots(len(samples) + len(short_read_samples) + draw_gene_track)
    axes = np.atleast_1d(axes)  # in case there was only one subplot

    if draw_gene_track:
        self.gene_track(ax=axes[0], x_range=x_range, draw_other_genes=draw_other_genes)

    for i, (sname, sidx) in enumerate(samples.items()):
        self.sashimi_plot(sidx, sname, axes[i + draw_gene_track], junctions_of_interest, x_range=x_range, **long_read_params)

    for i, (sname, sidx) in enumerate(short_read_samples.items()):
        self.sashimi_plot_short_reads(sidx, sname, axes[i + len(samples) + draw_gene_track], junctions_of_interest, x_range=x_range, **long_read_params)

    return f, axes


def sashimi_plot_short_reads(self, samples=None, title='short read coverage', ax=None, junctions_of_interest=None, x_range=None,
                             y_range=None, log_y=True,
                             jparams=None, min_cov_th=.001, high_cov_th=.05, text_width=.02, arc_type='both', text_height=1,
                             exon_color='green'):
    '''Draws short read Sashimi plot of the gene.

    The Sashimi plot depicts the genomic coverage from short read sequencing as blocks, and junction coverage as arcs.

    :param samples: Names of the short read samples to be depicted (as a list).
    :param title: Specify the title of the axis.
    :param ax: Specify the axis.
    :param junctions_of_interest: List of int pairs to define junctions of interest (which are highlighed in the plots)
    :param x_range: Genomic positions to specify the x range of the plot.
    :param y_range: Range for the coverage axis of the plot. Note to include space for the junction arcs.
        If not specified, the range will be determined automatically.
    :param log_y: Log scale for the coverage.
    :param jparams: Define the apperance of junctions, depending on their priority.
        A list with three dicts, defining parameters for low coverage junctions, high coverage junctions, and junctions of interest.
        For default values, see isotools._gene_plots.DEFAULT_JPARAMS
    :param exon_color: Specify the color of the genomic coverage blocks (e.g. the exons)
    :param high_cov_th: Minimum coverage for a junction to be considdered high coverage.
    :param min_cov_th: Coverage threshold for a junction to be considdered at all.
    :param text_width: Control the horizontal space that gets reserved for labels on the arcs. This affects the height of the arcs.
    :param arc_type: Label the junction arcs with  the "coverage" (e.g. number of supporting reads),
        "fraction" (e.g. fraction of supporting reads in %), or "both".
    :param text_height: Control the vertical space that gets reserved for labels on the arcs. This affects the height of the arcs.'''

    if samples is None:
        samples = list(self._transcriptome.infos['short_reads']['name'])  # all samples grouped # pylint: disable=W0212
    sidx = get_index(samples, self._transcriptome.infos['short_reads']['name'])  # pylint: disable=W0212
    if x_range is None:
        x_range = (self.start - 100, self.end + 100)

    if jparams is None:
        jparams = DEFAULT_JPARAMS

    short_reads = [self.short_reads(idx) for idx in sidx]
    # jparams=[low_cov_junctions,high_cov_junctions,interest_junctions]
    start = short_reads[0].reg[1]
    end = short_reads[0].reg[2]
    # delta=np.zeros(end-start)
    cov = np.zeros(end - start)
    junctions = {}
    for sr_cov in short_reads:
        cov += sr_cov.profile
        for k, v in sr_cov.junctions.items():
            junctions[k] = junctions.get(k, 0) + v
    if high_cov_th < 1:
        high_cov_th *= max(cov)
    if min_cov_th < 1:
        min_cov_th *= max(cov)
    if log_y:
        cov = np.log10(cov, where=cov > 0, out=np.nan * cov)
    # exons
    if ax is None:
        _, ax = plt.subplots()
    ax.fill_between(range(start, end), 0, cov, facecolor=exon_color)
    # junctions
    textpositions = []
    for (x1, x2), w in junctions.items():
        if junctions_of_interest is not None and (x1, x2) in junctions_of_interest:
            priority = 2
        elif w < min_cov_th:
            continue
        elif w < high_cov_th:
            priority = 0
        else:
            priority = 1
        y1 = cov[x1 - start - 1]
        y2 = cov[x2 - start]
        center = (x1 + x2) / 2
        width = x2 - x1
        bow_height = text_height
        if jparams[priority]['draw_label']:
            while any(_label_overlap((center, max(y1, y2) + bow_height), tp, text_width, text_height) for tp in textpositions):
                bow_height += text_height
            textpositions.append((center, max(y1, y2) + bow_height))
        if y1 < y2:
            bow_height = (y2 - y1 + bow_height, bow_height)
        elif y1 > y2:
            bow_height = (bow_height, bow_height + y1 - y2)
        else:
            bow_height = (bow_height, bow_height)
        bow1 = patches.Arc((center, y1), width=width, height=bow_height[0] * 2, theta1=90, theta2=180,
                           linewidth=jparams[priority]['lwd'], edgecolor=jparams[priority]['color'], zorder=priority)
        bow2 = patches.Arc((center, y2), width=width, height=bow_height[1] * 2, theta1=0, theta2=90,
                           linewidth=jparams[priority]['lwd'], edgecolor=jparams[priority]['color'], zorder=priority)
        ax.add_patch(bow1)
        ax.add_patch(bow2)
        if jparams[priority]['draw_label']:
            _ = ax.text(center, max(y1, y2) + min(bow_height) + text_height / 3, w, horizontalalignment='center', verticalalignment='bottom',
                        bbox=dict(boxstyle='round', facecolor='wheat', edgecolor=None, alpha=0.5)).set_clip_on(True)
        # bbox_list.append(txt.get_tightbbox(renderer = fig.canvas.renderer))

    ax.set_xlim(*x_range)
    if y_range is not None:
        ax.set_ylim(*y_range)
    if textpositions:
        ax.set_ylim(-text_height, max(tp[1] for tp in textpositions) + 2 * text_height)
    else:
        ax.set_ylim(-text_height, 3)  # todo: adjust y axis and ticklabels to coverage
    ax.set(frame_on=False)
    if log_y:
        ax.set_yticks([0, 1, 2, 3])
        ax.set_yticklabels([1, 10, 100, 1000])
    # ax.ticklabel_format(axis='x', style='sci',scilimits=(6,6))
    ax.set_title(title)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos=None: f'{x:,.0f}'))
    return ax


def sashimi_plot(self, samples=None, title='Long read sashimi plot', ax=None, junctions_of_interest=None, x_range=None, select_transcripts=None,
                 y_range=None, log_y=True,
                 jparams=None, exon_color='green', min_cov_th=.001, high_cov_th=.05, text_width=1,
                 arc_type='both', text_height=1):
    '''Draws long read Sashimi plot of the gene.

    The Sashimi plot depicts the genomic long read sequencing coverage of one or more samples as blocks, and junction coverage as arcs.


    :param samples: Names of the samples to be depicted (as a list).
    :param title: Specify the title of the axis.
    :param ax: Specify the axis.
    :param junctions_of_interest: List of int pairs to define junctions of interest (which are highlighed in the plots)
    :param x_range: Genomic positions to specify the x range of the plot.
        If not specified, the range will be determined to include the complete gene.
    :param y_range: Range for the coverage axis of the plot. Note to include space for the junction arcs.
        If not specified, the range will be determined automatically.
    :param log_y: Log scale for the coverage.
    :param select_transcripts: A list of transcript numbers from which the coverage is to be depicted.
        If obmitted, all transcripts are displayed.
    :param jparams: Define the apperance of junctions, depending on their priority.
        A list with three dicts, defining parameters for low coverage junctions, high coverage junctions, and junctions of interest.
        For default values, see isotools._gene_plots.DEFAULT_JPARAMS
    :param exon_color: Specify the color of the genomic coverage blocks (e.g. the exons)
    :param high_cov_th: Minimum coverage for a junction to be considdered high coverage.
    :param min_cov_th: Coverage threshold for a junction to be considdered at all.
    :param text_width: Scaling factor for the horizontal space that gets reserved for labels on the arcs.
        This affects the height of the arcs.
    :param arc_type: Label the junction arcs with  the "coverage" (e.g. number of supporting reads),
        "fraction" (e.g. fraction of supporting reads in %), or "both".
    :param text_height: Scaling factor for the vertical space that gets reserved for labels on the arcs.
        This affects the height of the arcs.'''

    sg = self.segment_graph
    if jparams is None:
        jparams = DEFAULT_JPARAMS
    if samples is None:
        samples = self._transcriptome.samples
    if ax is None:
        _, ax = plt.subplots()

    sidx = get_index(samples, self._transcriptome.samples)  # pylint: disable=W0212
    if junctions_of_interest is None:
        junctions_of_interest = []
    if x_range is None:
        x_range = self.start - 100, self.end + 100
    node_matrix = sg.get_node_matrix()
    if select_transcripts:
        try:
            _ = iter(select_transcripts)  # maybe only one transcript provided?
        except TypeError:
            select_transcripts = select_transcripts,
        mask = np.ones(node_matrix.shape[0], np.bool)
        mask[select_transcripts] = False
        node_matrix[mask, :] = 0
    boxes = [(node[0], node[1], self.coverage[np.ix_(sidx, node_matrix[:, i])].sum()) for i, node in enumerate(sg)]
    if log_y:
        boxes = [(s, e, log10(c) if c > 1 else c/10) for s, e, c in boxes]
    max_height = max(1, max(h for s, e, h in boxes if has_overlap(x_range, (s, e))))
    text_height = (max_height/10)*text_height
    text_width = (x_range[1] - x_range[0]) * .02 * text_width

    total_weight = self.coverage[sidx, :].sum()
    if high_cov_th < 1:
        high_cov_th = high_cov_th * total_weight
    if min_cov_th < 1:
        min_cov_th = min_cov_th * total_weight
    # idx=list(range(len(sg)))
    arcs = []
    for i, (_, ee, _, suc) in enumerate(sg):
        weights = {}
        for tr, next_i in suc.items():
            if sg[next_i][0] == ee or not has_overlap(x_range, (ee, sg[next_i][0])):
                continue
            if select_transcripts is not None and tr not in select_transcripts:
                continue
            tr_j_cov = self.coverage[np.ix_(sidx, [tr])].sum()
            if tr_j_cov:
                weights[next_i] = weights.get(next_i, 0)+tr_j_cov
        arcs_new = [(ee, boxes[i][2], sg[next_i][0], boxes[next_i][2], w) for next_i, w in weights.items()]
        if arcs_new:
            arcs.extend(arcs_new)
    if ax is None:
        _, ax = plt.subplots(1)

    for st, end, h in boxes:
        if h > 0 & has_overlap(x_range, (st, end)):
            rect = patches.Rectangle((st, 0), (end - st), h, linewidth=1, edgecolor=exon_color, facecolor=exon_color, zorder=5)
            ax.add_patch(rect)
    textpositions = []
    for x1, y1, x2, y2, w in arcs:
        if not has_overlap(x_range, (x1, x2)):
            continue
        if junctions_of_interest is not None and (x1, x2) in junctions_of_interest:
            priority = 2
        elif w < min_cov_th:
            continue
        elif w < high_cov_th:
            priority = 0
        else:
            priority = 1
        text_x = (x1 + x2) / 2
        textalign = 'center'
        if text_x > x_range[1]:
            text_x = x_range[1]
            textalign = 'right'
        elif text_x < x_range[0]:
            text_x = x_range[0]
            textalign = 'left'
        width = x2 - x1
        bow_height = text_height

        if jparams[priority]['draw_label']:
            while any(_label_overlap((text_x, max(y1, y2) + bow_height), tp, text_width, text_height) for tp in textpositions):
                bow_height += text_height
            textpositions.append((text_x, max(y1, y2) + bow_height))
        if y1 < y2:
            bow_height = (y2 - y1 + bow_height, bow_height)
        elif y1 > y2:
            bow_height = (bow_height, bow_height + y1 - y2)
        else:
            bow_height = (bow_height, bow_height)
        bow1 = patches.Arc(((x1 + x2) / 2, y1), width=width, height=bow_height[0] * 2, theta1=90, theta2=180,
                           linewidth=jparams[priority]['lwd'], edgecolor=jparams[priority]['color'], zorder=priority)
        bow2 = patches.Arc(((x1 + x2) / 2, y2), width=width, height=bow_height[1] * 2, theta1=0, theta2=90,
                           linewidth=jparams[priority]['lwd'], edgecolor=jparams[priority]['color'], zorder=priority)
        ax.add_patch(bow1)
        ax.add_patch(bow2)

        if jparams[priority]['draw_label']:
            if arc_type == 'coverage':
                lab = str(w)
            else:  # fraction
                lab = f'{w/total_weight:.1%}'
                if arc_type == 'both':
                    lab = str(w) + ' / ' + lab
            _ = ax.text(text_x, max(y1, y2) + min(bow_height) + text_height / 3, lab,
                        horizontalalignment=textalign, verticalalignment='bottom', zorder=10 + priority,
                        bbox=dict(boxstyle='round', facecolor='wheat', edgecolor=None, alpha=0.5)).set_clip_on(True)
        # bbox_list.append(txt.get_tightbbox(renderer = fig.canvas.renderer))
    if y_range:
        ax.set_ylim(*y_range)
    elif textpositions:
        ax.set_ylim(-text_height, max(tp[1] for tp in textpositions) + 2 * text_height)
    else:
        ax.set_ylim(-text_height, max_height+text_height)

    ax.set_xlim(*x_range)
    ax.set(frame_on=False)
    if log_y:
        ax.set_yticks([0, 1, 2, 3])
        ax.set_yticklabels([1, 10, 100, 1000])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos=None: f'{x:,.0f}'))
    # ax.ticklabel_format(axis='x', style='sci',scilimits=(6,6))
    # ax.set_xscale(1e-6, 'linear')
    ax.set_title(title)


def gene_track(self, ax=None, title=None, reference=True, select_transcripts=None, label_exon_numbers=True,
               label_transcripts=True, label_fontsize=10, color='blue', x_range=None, draw_other_genes=False,
               query=None, mincoverage=None, maxcoverage=None):
    '''Draws a gene track of the gene.

    The gene track depicts the exon structure of a gene, like in a genome browser.
    Exons are depicted as boxes, and junctions are lines. For coding regions, the height of the boxes is increased.
    Transcripts are labeled with the name, and a ">" or "<" sign, marking the direction of transcription.

    :param ax: Specify the axis.
    :param title: Specify the title of the axis.
    :param reference: If True, depict the reference transcripts, else transcripts are defined by long read sequencing.
    :param select_transcripts: A list of transcript numbers to be depicted.
        If draw_other_genes is set, select_transcripts should be a dict with gene_name as keys and lists with transcript numbers as values.
        If obmitted, all transcripts are displayed.
    :param label_exon_numbers: Draw exon numbers within exons.
    :param label_transcripts: Draw transcript name below transcripts.
    :param label_fontsize: Specify the font sice for the labels.
    :param color: Specify the color for the exons.
    :param x_range: Genomic positions to specify the x range of the plot.
    :param draw_other_genes: If set to True, transcripts from other genes overlapping the depicted region are also displayed.
        You can also provide a list of gene names/ids, to specify which other genes should be included.
    :param query: Filter query, which is passed to Gene.filter_transcripts or Gene.filter_ref_transcripts
    :param mincoverage: Minimum coverage for the transcript to be depeicted. Ignored in case of reference=True.
    :param maxcoverage: Maximum coverage for the transcript to be depeicted. Ignored in case of reference=True.'''

    if select_transcripts is None:
        select_transcripts = {}
    elif isinstance(select_transcripts, list):
        select_transcripts = {self.name: select_transcripts}
    else:
        try:
            _ = iter(select_transcripts)
        except TypeError:
            select_transcripts = {self.name: [select_transcripts]}
    contrast = 'white' if np.mean(plt_col.to_rgb(color)) < .5 else 'black'
    if ax is None:
        _, ax = plt.subplots(1)
    if x_range is None:
        x_range = (self.start - 100, self.end + 100)
    blocked = []
    if draw_other_genes:
        if isinstance(draw_other_genes, list):
            ol_genes = {self._transcriptome[g] for g in draw_other_genes}.add(self)
        else:
            ol_genes = self._transcriptome.data[self.chrom].overlap(*x_range)
    else:
        ol_genes = {self}
    transcript_list = []
    for g in ol_genes:
        select_tr = g.filter_ref_transcripts(query) if reference else g.filter_transcripts(query, mincoverage, maxcoverage)
        if select_transcripts.get(g.name):
            select_tr = [trid for trid in select_tr if trid in select_transcripts.get(g.name)]
        if reference:  # select transcripts and sort by start
            transcript_list.extend([(g, tr_nr, tr) for tr_nr, tr in enumerate(g.ref_transcripts) if tr_nr in select_tr])
        else:
            transcript_list.extend([(g, tr_nr, tr) for tr_nr, tr in enumerate(g.transcripts) if tr_nr in select_tr])
    transcript_list.sort(key=lambda x: x[2]['exons'][0][0])  # sort by start position
    for g, tr_nr, tr in transcript_list:
        tr_start, tr_end = tr['exons'][0][0], tr['exons'][-1][1]
        if (tr_end < x_range[0] or tr_start > x_range[1]):  # transcript does not overlap x_range
            continue
        trid = '> ' if g.strand == '+' else '< '  # indicate the strand like in ensembl browser
        trid += tr['transcript_name'] if 'transcript_name' in tr else f'{g.name}_{tr_nr}'

        # find next line that is not blocked
        try:
            i = next(idx for idx, last in enumerate(blocked) if last < tr['exons'][0][0])
        except StopIteration:
            i = len(blocked)
            blocked.append(tr_end)
        else:
            blocked[i] = tr_end
        # line from TSS to PAS at 0.25
        ax.plot((tr_start, tr_end), [i + .25] * 2, color=color)
        if label_transcripts:
            pos = (max(tr_start, x_range[0]) + min(tr_end, x_range[1])) / 2
            ax.text(pos, i - .02, trid, ha='center', va='top', fontsize=label_fontsize, clip_on=True)
        for j, (st, end) in enumerate(tr['exons']):
            cds = None
            if 'CDS' in tr or 'ORF' in tr:
                cds = tr['CDS'] if 'CDS' in tr else tr['ORF']
            if cds is not None and cds[0] <= end and cds[1] >= st:  # CODING exon
                c_st, c_end = max(st, cds[0]), min(cds[1], end)  # coding start and coding end
                if c_st > st:  # first noncoding part
                    rect = patches.Rectangle((st, i + .125), (c_st - st), .25, linewidth=1, edgecolor=color, facecolor=color)
                    ax.add_patch(rect)
                if c_end < end:  # 2nd noncoding part
                    rect = patches.Rectangle((c_end, i + .125), (end - c_end), .25, linewidth=1, edgecolor=color, facecolor=color)
                    ax.add_patch(rect)
                # Coding part
                rect = patches.Rectangle((c_st, i), (c_end - c_st), .5, linewidth=1, edgecolor=color, facecolor=color)
                ax.add_patch(rect)
            else:  # non coding
                rect = patches.Rectangle((st, i + .125), (end - st), .25, linewidth=1, edgecolor=color, facecolor=color)
                ax.add_patch(rect)
            if label_exon_numbers and (end > x_range[0] and st < x_range[1]):
                enr = j + 1 if g.strand == '+' else len(tr['exons']) - j
                pos = (max(st, x_range[0]) + min(end, x_range[1])) / 2
                ax.text(pos, i + .25, enr, ha='center', va='center', color=contrast, fontsize=label_fontsize,
                        clip_on=True)  # bbox=dict(boxstyle='round', facecolor='wheat',edgecolor=None,  alpha=0.5)
        i += 1
    if title is None:
        title = f'{self.name} ({self.region})'
    ax.set_title(title)
    ax.set(frame_on=False)
    ax.get_yaxis().set_visible(False)
    ax.set_ylim(-.5, len(blocked))
    ax.set_xlim(*x_range)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos=None: f'{x:,.0f}'))
    return ax


def find_blocks(pos, segments, remove_zero_gaps=False):
    adj_pos = []
    offset = 0
    idx = 0
    try:
        while pos[0] > offset+segments[idx][1] - segments[idx][0]:
            offset += segments[idx][1] - segments[idx][0]
            idx += 1
        adj_pos = [[segments[idx][0]+pos[0]-offset, None]]
        while pos[1] > offset+segments[idx][1] - segments[idx][0]:
            adj_pos[-1][1] = segments[idx][1]
            offset += segments[idx][1] - segments[idx][0]
            idx += 1
            adj_pos.append([segments[idx][0], None])
    except IndexError:
        logger.error(f'attempt to postitions {pos} blocks to segments {segments}')
        raise
    adj_pos[-1][1] = segments[idx][0]+pos[1]-offset
    if remove_zero_gaps:
        adj_pos_gaps = adj_pos
        adj_pos = []
        for seg in adj_pos_gaps:
            if not adj_pos or adj_pos[-1][1] < seg[0]:
                adj_pos.append(seg)
            else:
                adj_pos[-1][1] = seg[-1]
    return adj_pos


def get_rects(blocks, h=1, w=.1, connect=False,  **kwargs):
    rects = [patches.Rectangle((b[0], h), b[1]-b[0], w, **kwargs) for b in blocks]

    if connect:
        c = h+w/2
        rects.extend([patches.Polygon(np.array([[b1[1], c], [b2[0], c]]), closed=True, **kwargs) for b1, b2 in pairwise(blocks)])
    return (rects)


def find_segments(transcripts, orf_only=True, seperate_exons=False):
    '''Find exonic parts of the gene, with respect to trids.'''
    if orf_only:
        exon_list = []
        for tr in transcripts:
            cds_pos = tr.get('CDS', tr.get('ORF'))
            exon_list.append([])
            if cds_pos is None:
                continue
            for e in tr['exons']:
                if e[1] < cds_pos[0]:
                    continue
                if e[0] > cds_pos[1]:
                    break
                exon_list[-1].append([max(e[0], cds_pos[0]), min(e[1], cds_pos[1])])
    else:
        exon_list = [tr['exons'] for tr in transcripts]

    junctions = sorted([(pos, bool(j), i) for i, cds in enumerate(exon_list) for e in cds for j, pos in enumerate(e)])
    open_c = 0
    offset = 0
    genome_map = []  # genomic interval,  e.g.[([12345, 12445],0) ([12900-12980],100)]
    # can be used to calculate genomic position within blocks
    segments = [[] for _ in transcripts]
    pre_pos = None
    for pos, is_end, tr_i in junctions:
        if open_c > 0:
            offset += pos-pre_pos
        else:
            assert not is_end, f'more exons closed than opened before: {pos} at {junctions}'
            genome_map.append([pos, None])
        if not is_end:
            if seperate_exons or not segments[tr_i] or segments[tr_i][-1][1] < offset:
                segments[tr_i].append([offset, None])
            open_c += 1
        else:
            segments[tr_i][-1][1] = offset
            open_c -= 1
            if open_c == 0:
                genome_map[-1][1] = pos
        pre_pos = pos
    return segments, tuple(tuple(pos) for pos in genome_map)


def genome_pos_to_gene_segments(pos, genome_map, strict=True):
    pos = sorted(set(pos))
    offset = 0
    reverse_strand = genome_map[0][0] > genome_map[-1][1]
    if reverse_strand:
        genome_map = [(seg[1], seg[0]) for seg in reversed(genome_map)]
    mapped_pos = []
    i = 0
    for seg in genome_map:
        while seg[1] >= pos[i]:
            if seg[0] <= pos[i]:
                mapped_pos.append(offset+pos[i]-seg[0])
            elif not strict:
                mapped_pos.append(offset)
            else:
                mapped_pos.append(None)
            i += 1
            if i == len(pos):
                break
        else:
            offset += seg[1]-seg[0]
            continue
        break
    else:
        for i in range(i, len(pos)):
            mapped_pos.append(None if strict else offset)
    if reverse_strand:
        trlen = sum(seg[1]-seg[0] for seg in genome_map)
        mapped_pos = [trlen-mp if mp is not None else None for mp in mapped_pos]
    return {p: mp for p, mp in zip(pos, mapped_pos)}


def plot_domains(self, source, categories=None, trids=True, ref_trids=False, label='name', include_utr=False, seperate_exons=True,
                 x_ticks='gene', ax=None, domain_cols=DOMAIN_COLS,  max_overlap=5, highlight=None, highlight_col='red'):
    '''Plot exonic part of transcripts, together with protein domanis and annotations.

    :param source: Source of protein domains, e.g. "annotation", "hmmer" or "interpro", for domains added by the functions
        "add_annotation_domains", "add_hmmer_domains" or "add_interpro_domains" respectively.
    :param categories: List of domain categories to be depicted, default: all categories.
    :param trids: List of transcript indices to be depicted. If True/False, all/none transcripts are depicted.
    :param ref_trids: List of reference transcript indices to be depicted. If True/False, all/none reference transcripts are depicted.
    :param label: Specify the type of label: eiter None, or id, or name.
    :param include_utr: If set True, the untranslated regions are also depicted.
    :param seperate_exons: If set True, exon boundaries are marked.
    :param x_ticks: Either "gene" or "genome". If set to "gene", positions are relative to the gene (continuous, starting from 0).
        If set to "genome", positions are (discontinous) genomic coordinates.
    :param ax: Specify the axis.
    :param domain_cols: Dicionary for the colors of different domain types.
    :param max_overlap: Maximum number of overlapping domains to be depicted. Longer domains have priority over shorter domains.
    :param highlight: List of genomic positions or intervals to highlight.
    :param highlight_col: Specify the color for highlight positions.'''

    if label is not None:
        assert label in ('id', 'name'), 'label needs to be either "id" or "name" (or None).'
        label_idx = 0 if label == 'id' else 1
    domain_cols = {k.lower(): v for k, v in domain_cols.items()}
    assert x_ticks in ["gene", "genome"], f'x_ticks should be "gene" or "genome", not "{x_ticks}"'
    if isinstance(trids, bool):
        trids = list(range(len(self.transcripts))) if trids else []
    trids = [trid for trid in trids if 'ORF' in self.transcripts[trid] or 'CDS' in self.transcripts[trid]]
    if isinstance(ref_trids, bool):
        ref_trids = list(range(len(self.ref_transcripts))) if ref_trids else []
    ref_trids = [trid for trid in ref_trids if 'ORF' in self.ref_transcripts[trid] or 'CDS' in self.ref_transcripts[trid]]
    transcripts = [(i, self.ref_transcripts[i]) for i in ref_trids] + [(i, self.transcripts[i]) for i in trids]
    n_tr = len(transcripts)
    if not transcripts:
        logger.error('no transcripts with ORF specified')
        return
    if ax is None:
        _, ax = plt.subplots(1)
    skipped = 0
    segments, genome_map = find_segments([tr for _, tr in transcripts], orf_only=not include_utr, seperate_exons=seperate_exons)
    max_len = max(seg[-1][1] for seg in segments)
    assert max_len == sum(seg[1]-seg[0] for seg in genome_map)
    if self.strand == "-":
        segments = [[[max_len-pos[1], max_len-pos[0]] for pos in reversed(seg)] for seg in segments]
        genome_map = tuple((pos[1], pos[0]) for pos in reversed(genome_map))
    if highlight is not None:
        highlight_pos = set()
        for pos in highlight:
            if isinstance(pos, collections.abc.Sequence):
                highlight_pos.update(pos)
            else:
                highlight_pos.add(pos)
        highlight_pos = sorted(highlight_pos)
        pos_map = genome_pos_to_gene_segments(highlight_pos, genome_map, False)
        for pos in highlight:
            if isinstance(pos, collections.abc.Sequence):
                assert len(pos) == 2, 'provide intervals as a sequence of length 2'
                # draw box
                box_x = sorted(pos_map[p] for p in pos)
                patch = patches.Rectangle((box_x[0], -n_tr+.7), box_x[1]-box_x[0], n_tr+.6, edgecolor=highlight_col, facecolor=highlight_col)
                ax.add_patch(patch)
            else:  # draw line
                ax.vlines(pos_map[pos], -n_tr+.7, 0.8, colors=[highlight_col])

    for line, (trid, tr) in enumerate(transcripts):
        seg = segments[line]
        orf_pos = tr.get('CDS', tr.get('ORF'))[:2]
        orf = sorted(self.find_transcript_positions(trid, orf_pos, reference=line < len(ref_trids)))
        if include_utr:
            for rect in get_rects(find_blocks((0, orf[0]), seg),
                                  h=-line+.2, w=.1, connect=True, linewidth=1, edgecolor="black", facecolor="white"):
                ax.add_patch(rect)
            for rect in get_rects(find_blocks((orf[1], sum(e[1]-e[0] for e in tr['exons'])), seg),
                                  h=-line+.2, w=.1, connect=True, linewidth=1, edgecolor="black", facecolor="white"):
                ax.add_patch(rect)
        else:
            orf = (0, orf[1]-orf[0])
        orf_block = find_blocks((orf), seg)
        for rect in get_rects(orf_block, h=-line, w=.5, connect=True, linewidth=1, edgecolor="black", facecolor="white"):
            ax.add_patch(rect)

        domains = [dom for dom in tr.get('domain', {}).get(source, []) if categories is None or dom[2] in categories]
        # sort by length
        domains.sort(key=lambda x: x[3][1]-x[3][0], reverse=True)
        # get positions relative to segments
        dom_blocks = [find_blocks([p+orf[0] for p in dom[3]], seg, True) for dom in domains]
        dom_line = {}
        for idx, block in enumerate(dom_blocks):
            i = 0
            block_interval = (block[0][0], block[-1][1])
            while any(has_overlap(block_interval, b[1]) for b in dom_line.setdefault(i, [])):
                i += 1
                if i >= max_overlap:
                    skipped += 1
                    break
            else:
                dom_line[i].append((idx, block_interval))  # idx in length-sorted domains

        w = .4/max(len(dom_line), 1)
        for dom_l in dom_line:
            h = -line+w*(len(dom_line)//2 + (dom_l+1)//2 * (-1 if dom_l % 2 else 1))+.05
            for idx, bl in dom_line[dom_l]:
                dom = domains[idx]
                try:
                    for rect in get_rects(dom_blocks[idx], h=h, w=w, linewidth=1, edgecolor="black", facecolor=domain_cols.get(dom[2].lower(), "white")):
                        ax.add_patch(rect)
                except IndexError:
                    logger.error(f'cannot add patch for {dom_blocks[idx]}')
                    raise
                if label is not None:
                    ax.text((bl[0]+bl[1])/2, h+w/2, dom[label_idx], ha='center', va='center', color='black', clip_on=True)

    if skipped:
        logger.warning("skipped %s domains, consider increasing max_overlap parameter", skipped)
    ax.set_ylim(-len(transcripts)+.5, 1)
    ax.set_xlim(-10, max_len+10)
    if x_ticks == 'genome':
        xticks = [0]+list(np.cumsum([abs(seg[1]-seg[0]) for seg in genome_map]))
        xticklabels = [str(genome_map[0][0])]+[f'{seg[0][1]}|{seg[1][0]}' for seg in pairwise(genome_map)] + [str(genome_map[-1][1])]
        ax.set_xticks(ticks=xticks, labels=xticklabels)
    ax.set_yticks(ticks=[-i+.25 for i in range(len(transcripts))], labels=[tr.get('transcript_name', f'{self.name} {trid}') for trid, tr in transcripts])
    return ax, genome_map
