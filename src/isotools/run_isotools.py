import isotools
import matplotlib.pyplot as plt
import argparse
# import numpy as np
import pandas as pd
from isotools import Transcriptome
import isotools.plots
import logging
import sys


logger = logging.getLogger('run_isotools')


def main():
    parser = argparse.ArgumentParser(prog='isotools', description='process LRTS data with isotool')
    parser.add_argument('--anno', metavar='<file.gtf/gff/gff3[.gz]>', help='specify reference annotation')
    parser.add_argument('--genome', metavar='<file.fasta>', help='specify reference genome file')
    parser.add_argument('--samples', metavar='<samples.tsv>', help='add samples from sample tsv')
    parser.add_argument('--file_prefix', metavar='</output/directory/prefix>', default='./', help='Specify output path and prefix.')
    parser.add_argument('--file_suffix', metavar='<suffix>', help='Specify output sufix (not used for pickle).')
    parser.add_argument('--short_read_samples', metavar='<samples.csv>', help='Specify tsv with short read samples.')
    parser.add_argument('--force_recreate', help='reimport transcriptomes from alignments, even in presence of pickle file.', action='store_true')
    parser.add_argument('--no_pickle', help='Do not pickle the transcriptome for later use.', action='store_true')
    parser.add_argument('--progress_bar', help='Show the progress of individual tasks.', action='store_true')
    parser.add_argument("-l", "--log", dest="logLevel", default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL', None], help="Set the logging level.")
    parser.add_argument('--group_by', metavar='<column name>',
                        help='specify column used for grouping the samples. This applies to \
                            --qc_plots, --altsplice_stats, --diff, --diff_plots and --altsplice_plots',
                        default='name')
    parser.add_argument('--custom_filter_tag', metavar='<TAG="expression">', help='add custom filter tag', nargs='*')
    parser.add_argument('--filter_query', metavar='<"expression">', default='not (INTERNAL_PRIMING or RTTS)',
                        help='filter the transcripts used in gtf and table output')
    parser.add_argument('--qc_plots', help='make qc plots', action='store_true')
    parser.add_argument('--altsplice_stats', help='alternative splicing barplots', action='store_true')
    parser.add_argument('--transcript_table', help='make transcript_table', action='store_true')
    parser.add_argument('--gtf_out', help='make filtered gtf', action='store_true')
    parser.add_argument('--diff', metavar='<group1/group2>', nargs='*', help='perform differential splicing analysis')
    parser.add_argument('--diff_plots', metavar='<n>', type=int, help='make sashimi plots for <n> top differential genes')
    parser.add_argument('--altsplice_plots', metavar='<n>', type=int,
                        help='make sashimi plots for <n> top covered alternative spliced genes for each category')

    args = parser.parse_args()

    if args.logLevel:
        logging.basicConfig(level=getattr(logging, args.logLevel), format='%(asctime)s %(levelname)s: %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S')

    logger.debug(f'arguments: {args}')
    if args.file_suffix is None:
        file_suffix = ''
    else:
        file_suffix = '_'+args.file_suffix

    try:
        isoseq = load_isoseq(args)
    except ValueError as e:
        logger.error(e)
        parser.print_help(sys.stderr)
        exit(1)

    groups = isoseq.groups(args.group_by)
    logger.debug(f'sample group definition: {groups}')

    if args.short_read_samples:
        illu_samples = pd.read_csv(args.short_read_samples)
        isoseq.add_short_read_coverage(dict(zip(illu_samples['name'], illu_samples['file_name'])))

    illu_groups = {}
    if 'short_reads' in isoseq.infos:  # todo: make this optional/parameter --dont_use_short_reads
        for grp, sa in groups.items():
            if sa in isoseq.infos['short_reads']['name']:
                i = pd.Index(isoseq.infos['short_reads']['name']).get_loc(sa)
                illu_groups.setdefault(grp, []).append(i)
        logger.debug(f'illumina sample group definition: {illu_groups}\n{isoseq.infos["short_reads"]}')

    if args.custom_filter_tag is not None:
        for f_def in args.custom_filter_tag:
            tag, f_expr = f_def.split('=', 2)
            if tag not in isoseq.filter['transcript']:
                logger.info('adding new filter rule %s in transcript context', tag)
            isoseq.add_filter(tag, f_expr, context='transcript', update=True)

    if args.transcript_table:
        trtab_fn = f'{args.file_prefix}_transcripts{file_suffix}.csv'
        logger.info(f'writing transcript table to {trtab_fn}')
        df = isoseq.transcript_table(groups=groups, coverage=True, tpm=True, query=args.filter_query, progress_bar=args.progress_bar)
        df.to_csv(trtab_fn)

    if args.gtf_out:
        gtf_fn = f'{args.file_prefix}_transcripts{file_suffix}.gtf'
        isoseq.write_gtf(gtf_fn, use_gene_name=True, query=args.filter_query, progress_bar=args.progress_bar)

    if args.qc_plots:
        filter_plots(isoseq, groups, f'{args.file_prefix}_filter_stats{file_suffix}.png', args.progress_bar)
        transcript_plots(isoseq, groups,  f'{args.file_prefix}_transcript_stats{file_suffix}.png', args.progress_bar)

    if args.altsplice_stats:
        altsplice_plots(isoseq, groups, f'{args.file_prefix}_altsplice{file_suffix}.png', args.progress_bar)

    if args.altsplice_plots:
        examples = altsplice_examples(isoseq, args.altsplice_plots)
        plot_altsplice_examples(isoseq, groups,  illu_groups, examples, args.file_prefix, file_suffix)

    if args.diff is not None:
        test_differential(isoseq, groups, illu_groups, args, file_suffix)

    if not args.no_pickle:
        logger.info('saving transcripts as pickle file')
        isoseq.save(args.file_prefix+'_isotools.pkl')


def load_isoseq(args):
    isoseq = None
    # parameter checking:
    # if sample_tab is specified, genome must be specified
    if not args.force_recreate:
        try:
            isoseq = Transcriptome.load(args.file_prefix+'_isotools.pkl')
        except FileNotFoundError:
            if args.samples is None:
                raise ValueError('No samples specified')

    if args.samples:
        if args.anno is None or args.genome is None:
            raise ValueError('to add samples, genome and annotation must be provided.')
        if isoseq is None:
            isoseq = Transcriptome.from_reference(args.anno, progress_bar=args.progress_bar)
            isoseq.collapse_immune_genes()
        added = False
        sample_tab = pd.read_csv(args.samples, sep='\t')
        if 'sample_name' not in sample_tab.columns:
            logger.debug(sample_tab.columns)
            raise ValueError('No "sample_name" column found in sample table')
        if 'file_name' not in sample_tab.columns:
            raise ValueError('No "file_name" column found in sample table')
        for _, row in sample_tab.iterrows():
            if row['sample_name'] in isoseq.samples:
                logger.info('skipping already present sample %s', row["sample_name"])
                continue
            sample_args = {k: v for k, v in row.items() if k != 'file_name'}
            isoseq.add_sample_from_bam(fn=row.file_name, progress_bar=args.progress_bar, **sample_args)
            added = True
        if added:
            isoseq.add_qc_metrics(args.genome, progress_bar=args.progress_bar)
            isoseq.make_index()

    return isoseq


def filter_plots(isoseq, groups, filename, progress_bar):
    logger.info('filter statistics plots')
    f_stats = isoseq.filter_stats(groups=groups, weight_by_coverage=True, min_coverage=1, tr_filter={'progress_bar': progress_bar})
    plt.rcParams["figure.figsize"] = (15+5*len(groups), 7)
    fig, ax = plt.subplots()
    isotools.plots.plot_bar(f_stats[0], ax=ax, **f_stats[1])
    fig.tight_layout(rect=[0, 0, 1, .95])

    fig.savefig(filename)


def transcript_plots(isoseq, groups, filename, progress_bar):
    logger.info('perparing summary of quality control metrics...')
    logger.info('1) Number of RTTS, fragmentation and internal priming artefacts')
    f_stats = isoseq.filter_stats(groups=groups, weight_by_coverage=True, min_coverage=1,
                                  tr_filter={'progress_bar': progress_bar}, tags=('RTTS', 'FRAGMENT', 'INTERNAL_PRIMING'))
    tr_stats = []
    logger.info('2) Transcript length distribution')
    tr_stats.append(isoseq.transcript_length_hist(groups=groups, add_reference=True, min_coverage=2, tr_filter=dict(query='FSM', progress_bar=progress_bar)))
    logger.info('3) Distribution of downstream A fraction in known genes')
    tr_stats.append(isoseq.downstream_a_hist(groups=groups, tr_filter=dict(
        query='not (NOVEL_GENE or UNSPLICED)', progress_bar=progress_bar), ref_filter=dict(query='not UNSPLICED')))
    logger.info('4) Distribution of downstream A fraction in novel genes')
    tr_stats.append(isoseq.downstream_a_hist(groups=groups, tr_filter=dict(query='NOVEL_GENE and UNSPLICED', progress_bar=progress_bar)))
    logger.info('5) Distribution of direct repeats')
    tr_stats.append(isoseq.direct_repeat_hist(groups=groups, tr_filter=dict(progress_bar=progress_bar)))
    tr_stats.append((pd.concat([tr_stats[2][0].add_suffix(' novel unspliced'), tr_stats[1][0].add_suffix(' known multiexon')], axis=1), tr_stats[2][1]))

    plt.rcParams["figure.figsize"] = (30, 25)
    plt.rcParams.update({'font.size': 14})

    fig, axs = plt.subplots(3, 2)
    # A) transcript length
    isotools.plots.plot_distr(tr_stats[0][0], smooth=3, ax=axs[0, 0], **tr_stats[0][1])
    # D) frequency of artifacts
    isotools.plots.plot_bar(f_stats[0], ax=axs[0, 1], drop_categories=['PASS'], **f_stats[1])
    # B) internal priming
    isotools.plots.plot_distr(tr_stats[4][0][[c for c in tr_stats[4][0].columns if 'novel' in c]],
                              smooth=3, ax=axs[1, 0], density=True, fill=True, **tr_stats[4][1])
    isotools.plots.plot_distr(tr_stats[4][0][[c for c in tr_stats[4][0].columns if 'known' in c]],
                              smooth=3, ax=axs[1, 1], density=True, fill=True, **tr_stats[4][1])
    # C) RTTS
    isotools.plots.plot_distr(tr_stats[3][0][[c for c in tr_stats[3][0].columns if 'novel' in c]], ax=axs[2, 0], density=True, **tr_stats[3][1])
    isotools.plots.plot_distr(tr_stats[3][0][[c for c in tr_stats[3][0].columns if 'known' in c]], ax=axs[2, 1], density=True, **tr_stats[3][1])
    fig.tight_layout(rect=[0, 0, 1, .95])
    fig.savefig(filename)


def altsplice_plots(isoseq, groups, filename, progress_bar):
    logger.info('preparing novel splicing statistics...')
    altsplice = isoseq.altsplice_stats(groups=groups,  tr_filter=dict(query='not (RTTS or INTERNAL_PRIMING)', progress_bar=progress_bar))

    plt.rcParams["figure.figsize"] = (15+5*len(groups), 10)

    fig, ax = plt.subplots()
    isotools.plots.plot_bar(altsplice[0], ax=ax, drop_categories=['FSM'], **altsplice[1])
    fig.tight_layout(rect=[0, 0, 1, .95])
    fig.savefig(filename)


def altsplice_examples(isoseq, n, query='not FSM'):  # return the top n covered genes for each category
    examples = {}
    for g, trids, trs in isoseq.iter_transcripts(query=query, genewise=True):
        total_cov = g.coverage.sum()
        for trid, tr in zip(trids, trs):
            cov = g.coverage[:, trid].sum()
            score = cov*cov/total_cov
            for cat in tr['annotation'][1]:
                examples.setdefault(cat, []).append((score, g.name, g.id, trid, cov, total_cov))

    examples = {k: sorted(v, key=lambda x: -x[0]) for k, v in examples.items()}
    return {k: v[:n] for k, v in examples.items()}


def plot_altsplice_examples(isoseq, groups, illu_groups, examples, file_prefix, file_suffix):
    nplots = len(groups)+1
    # sample_idx = {r: i for i, r in enumerate(isoseq.infos['sample_table'].name)}
    if illu_groups:
        # illu_sample_idx = {r: i for i, r in enumerate(isoseq.infos['illumina_fn'])}
        if any(gn in illu_groups for gn in groups):
            illu_groups = {gn: illu_groups[gn] for gn in groups if gn in illu_groups}
        nplots += len(illu_groups)  # illumina is a dict with bam filenames

    plt.rcParams["figure.figsize"] = (20, 5*nplots)

    for cat, best_list in examples.items():
        logger.debug(cat+str(best_list))
        for i, (score, gene_name, gene_id, trid,  cov, total_cov) in enumerate(best_list):
            g = isoseq[gene_id]
            try:
                info = g.transcripts[trid]["annotation"][1][cat]
            except TypeError:
                info = list()
            logger.info(f'{i+1}. best example for {cat}: {gene_name} {trid} {info}, {cov} {total_cov} ({cov/total_cov:%})')
            joi = []  # set joi

            if info:
                junctions = []
                if cat == 'exon skipping':
                    exons = g.transcripts[trid]['exons']
                    for pos in info:
                        idx = next(i for i, e in enumerate(exons) if e[0] > pos[0])
                        junctions.append((exons[idx-1][1], exons[idx][0]))
                    info = junctions
                elif cat == 'novel exon':
                    exons = g.transcripts[trid]['exons']
                    for i, e in enumerate(exons[1:-1]):
                        if e in info:
                            junctions.extend([(exons[i][1], e[0]), (e[1], exons[i+2][0])])
                elif cat == 'novel junction':
                    junctions = info
                for pos in junctions:
                    try:
                        if len(pos) == 2 and all(isinstance(x, int) for x in pos):
                            joi.append(tuple(pos))  # if this is a junction, it gets highlighed in the plot
                    except TypeError:
                        pass
            print(f'junctions of interest: {joi}')
            fig, axs = g.sashimi_figure(samples=groups, junctions_of_interest=joi)

            fig.tight_layout()
            stem = f'{file_prefix}_altsplice{file_suffix}_{cat.replace(" ","_").replace("/","_")}_{g.name}'
            fig.savefig(f'{stem}_sashimi.png')
            # zoom
            if info:
                for pos in info:
                    if isinstance(pos, int):
                        start, end = pos, pos
                    elif len(pos) == 2 and all(isinstance(x, int) for x in pos):
                        if pos[1] < pos[0]:
                            start, end = sorted([pos[0], pos[0]+pos[1]])
                        else:
                            start, end = pos
                    else:
                        continue
                    for a in axs:
                        a.set_xlim((start - 100, end + 100))
                    axs[0].set_title(f'{g.name} {g.chrom}:{start}-{end} {cat} (cov={cov})')

                    plt.savefig(f'{stem}_zoom_{start}_{end}_sashimi.png')
            plt.close()


def plot_diffsplice(isoseq, de_tab, groups, illu_gr, file_prefix):

    nplots = len(groups)+1
    # sample_idx = {r: i for i, r in enumerate(isoseq.infos['sample_table'].name)}
    if illu_gr:
        # illu_sample_idx = {r: i for i, r in enumerate(isoseq.infos['illumina_fn'])}  # todo: add illumina
        nplots += len(illu_gr)

    plt.rcParams["figure.figsize"] = (20, 5*nplots)
    for gene_id in de_tab['gene_id'].unique():
        g = isoseq[gene_id]
        logger.info(f'sashimi plot for differentially spliced gene {g.name}')
        fig, axs = g.sashimi_figure(samples=groups)

        fig.savefig(f'{file_prefix}_{"_".join(groups)}_{g.name}_sashimi.png')
        # zoom
        for i, row in de_tab.loc[de_tab.gene == g.name].iterrows():
            if row.start > g.start and row.end < g.end:
                for a in axs:
                    a.set_xlim((row.start - 1000, row.end + 1000))
                axs[0].set_title(f'{g.name} {g.chrom}:{row.start}-{row.end}')
                fig.savefig(f'{file_prefix}_{"_".join(groups)}_{g.name}_zoom_{row.start}_{row.end}_sashimi.png')
        plt.close(fig)


def test_differential(isoseq, groups, illu_groups, args, file_suffix):
    file_prefix = f'{args.file_prefix}_diff{file_suffix}'
    for diff_cmp in args.diff:
        gr = diff_cmp.split('/')
        logger.debug(f'processing {gr}')
        if len(gr) != 2:
            logger.error('--diff argument format error: provide two groups seperated by "/" -- skipping')
            continue
        if not all(gn in groups for gn in gr):
            logger.error(
                f'--diff argument format error: group names {[gn for gn in gr if gn not in groups]} not found in sample table -- skipping')
            continue
        gr = {gn: groups[gn] for gn in gr}
        res = isoseq.altsplice_test(gr, progress_bar=args.progress_bar).sort_values('pvalue')
        sig = res.padj < .1
        logger.info(f'{sum(sig)} differential splice sites in {len(res.loc[sig,"gene"].unique())} genes for {" vs ".join(gr)}')
        res.to_csv(f'{file_prefix}_{"_".join(gr)}.csv', index=False)
        if args.diff_plots is not None:

            sig_tab = res.head(args.diff_plots)
            if illu_groups:
                illu_cov = list()
                for g, jstart, jend in zip(sig_tab.gene, sig_tab.start, sig_tab.end):
                    ji = (jstart, jend)
                    # g_cov=[0,0,0,0]
                    j_cov = [{}, {}]
                    cov = isoseq[g].illumina_coverage
                    for gi, grp_n in enumerate(gr):
                        if grp_n not in illu_groups:
                            j_cov[gi] = 'NA'
                        for sn in illu_groups[grp_n]:
                            i = illu_groups[sn]
                            for k, v in cov[i].junctions.items():
                                j_cov[gi][k] = j_cov[gi].get(k, 0)+v
                            j_cov[gi].setdefault(ji, 0)
                    illu_cov.append((j_cov[0][ji], j_cov[1][ji], max(j_cov[0].values()), max(j_cov[1].values())))
                illu_cov = {k: v for k, v in zip(['illu_cov1', 'illu_cov2', 'illu_max1', 'illu_max2'], zip(*illu_cov))}
                sig_tab = sig_tab.assign(**illu_cov)
            sig_tab.to_csv(f'{file_prefix}_top_{"_".join(gr)}.csv')
            plot_diffsplice(isoseq, res.head(args.diff_plots), gr, illu_groups, file_prefix)


if __name__ == '__main__':
    main()
