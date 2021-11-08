
import matplotlib.colors as plt_col
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
from math import log10
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


def sashimi_figure(self, samples=None, short_read_samples=None, draw_gene_track=True,
                   long_read_params=None, short_read_params=None, junctions_of_interest=None, x_range=None):
    '''Arranges multiple Sashimi plots of the gene.

    The Sashimi figure consist of a reference gene track, long read sashimi plots for one or more samples or groups of samples,
    and optionally short read sashimi plots for one or more samples or groups of samples.

    :param samples: Definition of samples (as a list) or groups of samples (as a dict) for long read plots.
    :param short_read_samples: Definition of samples (as a list) or groups of samples (as a dict) for short read plots.
    :param draw_gene_track: Specify whether to plot the reference gene track.
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
        self.gene_track(ax=axes[0], x_range=x_range)

    for i, (sname, sidx) in enumerate(samples.items()):
        self.sashimi_plot(sidx, sname, axes[i + draw_gene_track], junctions_of_interest, x_range=x_range, **long_read_params)

    for i, (sname, sidx) in enumerate(short_read_samples.items()):
        self.sashimi_plot_short_reads(sidx, sname, axes[i + len(samples) + draw_gene_track], junctions_of_interest, x_range=x_range, **long_read_params)

    return f, axes


def sashimi_plot_short_reads(self, samples=None, title='short read coverage', ax=None, junctions_of_interest=None, x_range=None,
                             jparams=None, min_cov_th=.001, high_cov_th=.05, text_width=.02, arc_type='both', text_height=1,
                             exon_color='green'):
    '''Draws short read Sashimi plot of the gene.

    The Sashimi plot depicts the genomic coverage from short read sequencing as blocks, and junction coverage as arcs.

    :param samples: Names of the short read samples to be depicted (as a list).
    :param title: Specify the title of the axis.
    :param ax: Specify the axis.
    :param junctions_of_interest: List of int pairs to define junctions of interest (which are highlighed in the plots)
    :param x_range: Genomic positions to specify the x range of the plot.
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
    # exons
    if ax is None:
        _, ax = plt.subplots()
    ax.fill_between(range(start, end), 0, np.log10(cov, where=cov > 0, out=np.nan * cov), facecolor=exon_color)
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
        y1 = cov[x1 - start - 1] + .5
        y2 = cov[x2 - start] + .5
        center = (x1 + x2) / 2
        width = x2 - x1
        bow_height = text_height
        while any(_label_overlap((center, log10(max(y1, y2)) + bow_height), tp, text_width, text_height) for tp in textpositions):
            bow_height += text_height
        textpositions.append((center, log10(max(y1, y2)) + bow_height))
        if y1 < y2:
            bow_height = (log10(y2 / y1) + bow_height, bow_height)
        elif y1 > y2:
            bow_height = (bow_height, bow_height + log10(y1 / y2))
        else:
            bow_height = (bow_height, bow_height)
        bow1 = patches.Arc((center, log10(y1)), width=width, height=bow_height[0] * 2, theta1=90, theta2=180,
                           linewidth=jparams[priority]['lwd'], edgecolor=jparams[priority]['color'], zorder=priority)
        bow2 = patches.Arc((center, log10(y2)), width=width, height=bow_height[1] * 2, theta1=0, theta2=90,
                           linewidth=jparams[priority]['lwd'], edgecolor=jparams[priority]['color'], zorder=priority)
        ax.add_patch(bow1)
        ax.add_patch(bow2)
        if jparams[priority]['draw_label']:
            _ = ax.text(center, log10(max(y1, y2)) + min(bow_height) + text_height / 3, w, horizontalalignment='center', verticalalignment='bottom',
                        bbox=dict(boxstyle='round', facecolor='wheat', edgecolor=None, alpha=0.5)).set_clip_on(True)
        # bbox_list.append(txt.get_tightbbox(renderer = fig.canvas.renderer))

    ax.set_xlim(*x_range)
    if textpositions:
        ax.set_ylim(-text_height, max(tp[1] for tp in textpositions) + 2 * text_height)
    else:
        ax.set_ylim(-text_height, 3)  # todo: adjust y axis and ticklabels to coverage
    ax.set(frame_on=False)
    ax.set_yticks([0, 1, 2, 3])
    ax.set_yticklabels([1, 10, 100, 1000])
    # ax.ticklabel_format(axis='x', style='sci',scilimits=(6,6))
    ax.set_title(title)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos=None: f'{x:,.0f}'))
    return ax


def sashimi_plot(self, samples=None, title='Long read sashimi plot', ax=None, junctions_of_interest=None, x_range=None, select_transcripts=None,
                 jparams=None, exon_color='green', min_cov_th=.001, high_cov_th=.05, text_width=.02,
                 arc_type='both', text_height=1):
    '''Draws long read Sashimi plot of the gene.

    The Sashimi plot depicts the genomic long read sequencing coverage of one or more samples as blocks, and junction coverage as arcs.


    :param samples: Names of the samples to be depicted (as a list).
    :param title: Specify the title of the axis.
    :param ax: Specify the axis.
    :param junctions_of_interest: List of int pairs to define junctions of interest (which are highlighed in the plots)
    :param x_range: Genomic positions to specify the x range of the plot.
    :param select_transcripts: A list of transcript numbers from which the coverage is to be depicted.
        If obmitted, all transcripts are displayed.
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

    if text_width < 1:
        text_width = (sg[-1][1] - sg[0][0]) * text_width
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
            if select_transcripts is not None and tr not in select_transcripts:
                continue
            weights.setdefault(next_i, 0)
            weights[next_i] += self.coverage[np.ix_(sidx, [tr])].sum()
        arcs_new = [(ee, boxes[i][2], sg[next_i][0], boxes[next_i][2], w) for next_i, w in weights.items() if sg[next_i][0] > ee and w > 0]
        if arcs_new:
            arcs.extend(arcs_new)
    if ax is None:
        _, ax = plt.subplots(1)

    for st, end, h in boxes:
        if h > 0:
            rect = patches.Rectangle((st, 0), (end - st), log10(h), linewidth=1, edgecolor=exon_color, facecolor=exon_color, zorder=5)
            ax.add_patch(rect)
    textpositions = []
    for x1, y1, x2, y2, w in arcs:
        if junctions_of_interest is not None and (x1, x2) in junctions_of_interest:
            priority = 2
        elif w < min_cov_th:
            continue
        elif w < high_cov_th:
            priority = 0
        else:
            priority = 1
        center = (x1 + x2) / 2
        width = x2 - x1
        bow_height = text_height
        while any(_label_overlap((center, log10(max(y1, y2)) + bow_height), tp, text_width, text_height) for tp in textpositions):
            bow_height += text_height
        textpositions.append((center, log10(max(y1, y2)) + bow_height))
        if y1 < y2:
            bow_height = (log10(y2 / y1) + bow_height, bow_height)
        elif y1 > y2:
            bow_height = (bow_height, bow_height + log10(y1 / y2))
        else:
            bow_height = (bow_height, bow_height)
        bow1 = patches.Arc((center, log10(y1)), width=width, height=bow_height[0] * 2, theta1=90, theta2=180,
                           linewidth=jparams[priority]['lwd'], edgecolor=jparams[priority]['color'], zorder=priority)
        bow2 = patches.Arc((center, log10(y2)), width=width, height=bow_height[1] * 2, theta1=0, theta2=90,
                           linewidth=jparams[priority]['lwd'], edgecolor=jparams[priority]['color'], zorder=priority)
        ax.add_patch(bow1)
        ax.add_patch(bow2)
        if arc_type == 'coverage':
            lab = str(w)
        else:  # fraction
            lab = f'{w/total_weight:.1%}'
            if arc_type == 'both':
                lab = str(w) + ' / ' + lab
        if jparams[priority]['draw_label']:
            _ = ax.text(center, log10(max(y1, y2)) + min(bow_height) + text_height / 3, lab,
                        horizontalalignment='center', verticalalignment='bottom', zorder=10 + priority,
                        bbox=dict(boxstyle='round', facecolor='wheat', edgecolor=None, alpha=0.5)).set_clip_on(True)
        # bbox_list.append(txt.get_tightbbox(renderer = fig.canvas.renderer))
    if textpositions:
        ax.set_ylim(-text_height, max(tp[1] for tp in textpositions) + 2 * text_height)
    else:
        ax.set_ylim(-text_height, 3)  # todo: adjust y axis and ticklabels to coverage
    ax.set_xlim(*x_range)
    ax.set(frame_on=False)
    ax.set_yticks([0, 1, 2, 3])
    ax.set_yticklabels([1, 10, 100, 1000])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos=None: f'{x:,.0f}'))
    # ax.ticklabel_format(axis='x', style='sci',scilimits=(6,6))
    # ax.set_xscale(1e-6, 'linear')
    ax.set_title(title)


def gene_track(self, ax=None, title=None, reference=True, select_transcripts=None, label_exon_numbers=True,
               label_transcripts=True, label_fontsize=10, color='blue', x_range=None, draw_other_genes=False):
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
        You can also provide a list of gene names/ids, to specify which other genes should be included.'''

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
        select_tr = select_transcripts.get(g.name, None)
        if reference:  # select transcripts and sort by start
            transcript_list.extend([(g, tr_nr, tr) for tr_nr, tr in enumerate(g.ref_transcripts) if select_tr is None or tr_nr in select_tr])
        else:
            transcript_list.extend([(g, tr_nr, tr) for tr_nr, tr in enumerate(g.transcripts) if select_tr is None or tr_nr in select_tr])
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
            if 'CDS' in tr and tr['CDS'][0] <= end and tr['CDS'][1] >= st:  # CODING exon
                c_st, c_end = max(st, tr['CDS'][0]), min(tr['CDS'][1], end)  # coding start and coding end
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
