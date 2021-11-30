from scipy.stats import binom, norm, chi2, betabinom  # pylint: disable-msg=E0611
from scipy.special import gammaln, polygamma  # pylint: disable-msg=E0611
from scipy.optimize import minimize
import statsmodels.stats.multitest as multi
import logging
import numpy as np
import pandas as pd
import itertools

# from ._utils import overlap
# from .decorators import deprecated, debug, experimental
from ._utils import _filter_function

logger = logging.getLogger('isotools')

# differential splicing


def proportion_test(x, n):
    # Normal approximation
    # x,n should be lenght 2(the two groups)
    # tests H0: proportions are equal vs H1: proportions are different (two sided)
    x = [xi.sum() for xi in x]
    n = [ni.sum() for ni in n]
    p1 = [x[i] / n[i] for i in range(2)]
    p0 = (x[0] + x[1]) / (n[0] + n[1])
    z = abs(p1[0] - p1[1]) / np.sqrt(p0 * (1 - p0) * (1 / n[0] + 1 / n[1]))
    return(2 * norm.sf(z)), (p1[0], 0, p1[1], 0, p0, 0)  # two sided alternative


def binom_lr_test(x, n):
    # likelihood ratio test
    # x,n should be length 2 (the two groups)
    # principle: log likelihood ratio of M0/M1 is chi2 distributed
    x = [xi.sum() for xi in x]
    n = [ni.sum() for ni in n]
    p1 = [x[i] / n[i] for i in range(2)]
    p0 = (x[0] + x[1]) / (n[0] + n[1])
    # calculate the log likelihoods
    l0 = binom.logpmf(x, n, p0).sum()
    l1 = binom.logpmf(x, n, p1).sum()
    # calculate the pvalue (sf=1-csf(), 1df)
    return chi2.sf(2 * (l1 - l0), 1), (p1[0], 0, p1[1], 0, p0, 0)


def loglike_betabinom(params, k, n):
    '''returns  log likelihood of betabinomial and its partial derivatives'''
    a, b = params
    logpdf = gammaln(n + 1) + gammaln(k + a) + gammaln(n - k + b) + gammaln(a + b) - \
        (gammaln(k + 1) + gammaln(n - k + 1) + gammaln(a) + gammaln(b) + gammaln(n + a + b))
    e = polygamma(0, a + b) - polygamma(0, n + a + b)
    da = e + polygamma(0, k + a) - polygamma(0, a)
    db = e + polygamma(0, n - k + b) - polygamma(0, b)
    return -np.sum(logpdf), np.array((-np.sum(da), -np.sum(db)))


def betabinom_ml(xi, ni):
    '''Calculate maximum likelihood parameter of beta binomial distribution for a group of samples with xi successes and ni trials.

    :param xi: number of successes, here coverage of the alternative for all samples of the group as 1d numpy array
    :param ni: number of trials, here total coverage for the two sample groups for all samples of the group as 1d numpy array
    '''
    # x and n must be np arrays
    if sum(ni) == 0:
        params = params_alt = None, None
        return params, params_alt, False
    xi, ni = xi[ni > 0], ni[ni > 0]  # avoid div by 0
    prob = xi / ni
    m = prob.mean()  # estimate initial parameters
    d = prob.var()
    success = True
    if d == 0:  # just one sample? or all exactly the same proportion
        params = params_alt = m, None  # in this case the betabinomial reduces to the binomial
    else:
        d = max(d, 1e-6)  # to avoid division by 0
        e = (m**2 - m + d)  # helper
        # find ml estimates for a and b
        mle = minimize(loglike_betabinom, x0=[-m * e / d, ((m - 1) * e) / d], bounds=((1e-6, None), (1e-6, None)),
                       args=(xi, ni), options={'maxiter': 250}, method='L-BFGS-B', jac=True)
        a, b = params = mle.x
        params_alt = (a / (a + b), a * b / ((a + b)**2 * (a + b + 1)))  # get alternative parametrization (mu and disp)
        # mle = minimize(loglike_betabinom2, x0=[-d/(m*e),d/((m-1)*e)],bounds=((1e-9,None),(1e-9,None)),
        #   args=(xi,ni),options={'maxiter': 250}, method='L-BFGS-B', tol=1e-6)
        # params=([1/p for p in mle.x])
        params_alt = (a / (a + b), a * b / ((a + b)**2 * (a + b + 1)))  # get alternative parametrization (mu and disp)

        if not mle.success:
            # should not happen to often, mainly with mu close to boundaries
            logger.debug(f'no convergence in betabinomial fit: k={xi}\nn={ni}\nparams={params}\nmessage={mle.message}')
            success = False  # prevent calculation of p-values based on non optimal parameters
    return params, params_alt, success


def betabinom_lr_test(x, n):
    ''' Likelihood ratio test with random-effects betabinomial model.

    This test modles x as betabinomial(n,a,b), eg a binomial distribution, where p follows beta ditribution with parameters a,b>0
    mean m=a/(a+b) overdispersion d=ab/((a+b+1)(a+b)^2) --> a=-m(m^2-m+d)/d b=(m-1)(m^2-m+d)/d
    principle: log likelihood ratio of M0/M1 is chi2 distributed

    :param x: coverage of the alternative for the two sample groups
    :param n: total coverage for the two sample groups'''

    if any(ni.sum() == 0 for ni in n):
        return (np.nan, [None, None])  # one group is not covered at all - no test possible. Checking this to avoid RuntimeWarnings (Mean of empty slice)
    x_all, n_all = (np.concatenate(x), np.concatenate(n))
    # calculate ml parameters
    ml_1 = betabinom_ml(x[0], n[0])
    ml_2 = betabinom_ml(x[1], n[1])
    ml_all = betabinom_ml(x_all, n_all)

    if not (ml_1[2] and ml_2[2] and ml_all[2]):  # check success
        return np.nan, list(ml_1[1] + ml_2[1] + ml_all[1])
    try:
        l0 = betabinom_ll(x_all, n_all, *ml_all[0]).sum()
        l1 = betabinom_ll(x[0], n[0], *ml_1[0]).sum() + betabinom_ll(x[1], n[1], *ml_2[0]).sum()
    except (ValueError, TypeError):
        logger.critical(f'betabinom error: x={x}\nn={n}\nparams={ml_1[0]}/{ml_2[0]}/{ml_all[0]}')  # should not happen
        raise
    return chi2.sf(2 * (l1 - l0), 2), list(ml_1[1] + ml_2[1] + ml_all[1])  # note that we need two degrees of freedom here as h0 hsa two parameters, h1 has 4


def betabinom_ll(x, n, a, b):
    if b is None:
        return binom.logpmf(x, n, a).sum()
    else:
        return betabinom.logpmf(x, n, a, b).sum()


TESTS = {'betabinom_lr': betabinom_lr_test,
         'binom_lr': binom_lr_test,
         'proportions': proportion_test}


def altsplice_test(self, groups, min_total=100, min_alt_fraction=.1, min_n=10, min_sa=.51, test='auto', padj_method='fdr_bh', types=None, progress_bar=True):
    '''Performs the alternative splicing event test.

    :param groups: Dict with groupnames as keys and lists of samplenames as values, defining the two groups for the test.
        If more then two groups are provided, test is performed between first two groups, but maximum likelihood parameters
        (expected PSI and dispersion) will be computet for the other groups as well.
    :param min_total: Minimum total coverage over all selected samples (for both groups combined).
    :param min_alt_fraction: Minimum fraction of reads supporting the alternative (for both groups combined).
    :param min_n: The minimum coverage of the event for an individual sample to be considered for the min_sa filter.
    :param min_sa: The fraction of samples within each group that must be covered by at least min_n reads.
    :param test: The name of one of the implemented statistical tests ('betabinom_lr','binom_lr','proportions').
    :param padj_method: Specify the method for multiple testing correction.
    :param types: Restrict the analysis on types of events. If ommited, all types are tested.'''
    # assert len(groups) == 2 , "length of groups should be 2, but found %i" % len(groups)
    # find groups and sample indices
    if isinstance(groups, dict):
        groupnames = list(groups)
        groups = list(groups.values())
    elif all(isinstance(gn, str) and gn in self.groups() for gn in groups):
        groupnames = list(groups)
        groups = [self.groups()[gn] for gn in groupnames]
    elif all(isinstance(grp, list) for grp in groups):
        groupnames = ['group{i+1}' for i in range(len(groups))]
    else:
        raise ValueError('groups not found in dataset')
    notfound = [sa for grp in groups for sa in grp if sa not in self.samples]
    if notfound:
        raise ValueError(f"Cannot find the following samples: {notfound}")
    assert not (groupnames[0] in groupnames[1] or groupnames[1] in groupnames[0]), 'group names must not be contained in other group names'
    if isinstance(test, str):
        if test == 'auto':
            test = 'betabinom_lr' if min(len(g) for g in groups[:2]) > 1 else 'proportions'
        test_name = test
        try:
            test = TESTS[test]
        except KeyError as e:
            raise ValueError(f'test must be one of {str(list(TESTS))}') from e
    else:
        test_name = 'custom'

    logger.info('testing differential splicing for %s using %s test', ' vs '.join(f'{groupnames[i]} ({len(groups[i])})' for i in range(2)), test_name)
    sa_idx = {sa: idx[0] for sa, idx in self._get_sample_idx().items()}
    grp_idx = [[sa_idx[sa] for sa in grp] for grp in groups]
    sidx = grp_idx[0] + grp_idx[1]
    if min_sa < 1:
        min_sa *= sum(len(gr) for gr in groups[:2])
    res = []
    for g in self.iter_genes(progress_bar=progress_bar):
        if g.coverage[sidx, :].sum() < min_total:
            continue
        known = {}  # check for known events
        if g.is_annotated and g.n_transcripts:
            sg = g.ref_segment_graph
            # find annotated alternatives for gene (e.g. known events)
            for _, _, nX, nY, splice_type in sg.find_splice_bubbles(types=types):
                if splice_type in ("TSS", "PAS"):
                    if (splice_type == "TSS") == (g.strand == "+"):
                        known.setdefault(splice_type, set()).add((sg[nX].end))
                    else:
                        known.setdefault(splice_type, set()).add((sg[nY].start))
                else:
                    known.setdefault(splice_type, set()).add((sg[nX].end, sg[nY].start))
        sg = g.segment_graph
        for setA, setB, nX, nY, splice_type in sg.find_splice_bubbles(types=types):

            junction_cov = g.coverage[:, setB].sum(1)
            total_cov = g.coverage[:, setA].sum(1) + junction_cov
            if total_cov[sidx].sum() < min_total:
                continue
            alt_fraction = junction_cov[sidx].sum() / total_cov[sidx].sum()
            if alt_fraction < min_alt_fraction or alt_fraction > 1 - min_alt_fraction:
                continue
            x = [junction_cov[grp] for grp in grp_idx]
            n = [total_cov[grp] for grp in grp_idx]
            if sum((ni >= min_n).sum() for ni in n[:2]) < min_sa:
                continue
            pval, params = test(x[:2], n[:2])
            params_other = tuple(v for xi, ni in zip(x[2:], n[2:]) for v in betabinom_ml(xi, ni)[1])
            if splice_type in ['TSS', 'PAS']:
                start, end = sg[nX].start, sg[nY].end
                if (splice_type == "TSS") == (g.strand == "+"):
                    novel = end not in known.get(splice_type, set())
                else:
                    novel = start not in known.get(splice_type, set())
            else:
                start, end = sg[nX].end, sg[nY].start
                novel = (start, end) not in known.get(splice_type, set())
            res.append(tuple(itertools.chain((g.name, g.id, g.chrom, g.strand, start, end, splice_type, novel, pval, setA, setB), params, params_other,
                                             (val for lists in zip(x, n) for pair in zip(*lists) for val in pair))))

    df = pd.DataFrame(res, columns=(['gene', 'gene_id', 'chrom', 'strand', 'start', 'end', 'splice_type', 'novel', 'pvalue', 'trA', 'trB'] +
                                    [gn + part for gn in groupnames[:2] + ['total'] + groupnames[2:] for part in ['_PSI', '_disp']] +
                                    [f'{sa}_{gn}_{w}' for gn, grp in zip(groupnames, groups) for sa in grp for w in ['in_cov', 'total_cov']]))
    try:
        mask = np.isfinite(df['pvalue'])
        padj = np.empty(mask.shape)
        padj.fill(np.nan)
        padj[mask] = multi.multipletests(df.loc[mask, 'pvalue'], method=padj_method)[1]
        df.insert(8, 'padj', padj)
    except TypeError as e:  # apparently this happens if df is empty...
        logger.error(f'unexpected error during calculation of adjusted p-values: {e}')
    return df


def alternative_splicing_events(self, min_total=100, min_alt_fraction=.1, samples=None, region=None, query=None, progress_bar=False):
    '''Finds alternative splicing events.

    Finds alternative splicing events and potential transcription start sites/polyA sites
    by searching for splice bubbles in the Segment Graph.
    Genes may be specified by genomic "region", and/or by filter tags / novelity class using the "query" parameters.

    :param min_total: Minimum total coverage over all selected samples.
    :param min_alt_fraction: Minimum fraction of reads supporting the alternative.
    :param samples: Specify the samples to consider. If omitted, all samples are selected.
    :param region: Specify the region, either as (chr, start, end) tuple or as "chr:start-end" string.
        If omitted, the complete genome is searched.
    :param query: Specify gene filter query.
    :param progress_bar: If set True the progress is shown.
    :return: Table with alternative splicing events.'''
    bubbles = []
    if samples is None:
        samples = self.samples
    assert all(s in self.samples for s in samples), 'not all specified samples found'
    sa_dict = {sa: i for i, sa in enumerate(self.samples)}
    sidx = np.array([sa_dict[sa] for sa in samples])

    assert 0 < min_alt_fraction < .5, 'min_alt_fraction must be > 0 and < 0.5'
    for g in self.iter_genes(region=region, query=query, progress_bar=progress_bar):
        if g.coverage[sidx, :].sum() < min_total:
            continue
        known = {}  # check for known events
        if g.is_annotated and g.n_transcripts:
            sg = g.ref_segment_graph
            for _, _, nX, nY, splice_type in sg.find_splice_bubbles():  # find annotated alternatives (known)
                if splice_type in ("TSS", "PAS"):
                    if (splice_type == "TSS") == (g.strand == "+"):
                        known.setdefault(splice_type, set()).add((sg[nX].end))
                    else:
                        known.setdefault(splice_type, set()).add((sg[nY].start))
                else:
                    known.setdefault(splice_type, set()).add((sg[nX].end, sg[nY].start))
        sg = g.segment_graph
        for setA, setB, nX, nY, splice_type in sg.find_splice_bubbles():
            junction_cov = g.coverage[np.ix_(sidx, setA)].sum(1)
            total_cov = g.coverage[np.ix_(sidx, setB)].sum(1) + junction_cov
            if total_cov.sum() >= min_total and min_alt_fraction < junction_cov.sum() / total_cov.sum() < 1 - min_alt_fraction:
                if splice_type in ['TSS', 'PAS']:
                    start, end = sg[nX].start, sg[nY].end
                    if (splice_type == "TSS") == (g.strand == "+"):
                        novel = end not in known.get(splice_type, set())
                    else:
                        novel = start not in known.get(splice_type, set())
                else:
                    start, end = sg[nX].end, sg[nY].start
                    novel = (start, end) not in known.get(splice_type, set())
                bubbles.append([g.id, g.chrom, start, end, splice_type, novel] + list(junction_cov) + list(total_cov))
    return pd.DataFrame(bubbles, columns=['gene', 'chr', 'start', 'end', 'splice_type', 'novel'] +
                        [f'{sa}_{what}' for what in ['in_cov', 'total_cov'] for sa in samples])

# summary tables (can be used as input to plot_bar / plot_dist)


def altsplice_stats(self, groups=None, weight_by_coverage=True, min_coverage=2, tr_filter={}):
    '''Summary statistics for novel alternative splicing.

    This function counts the novel alternative splicing events of LRTS isoforms with respect to the reference annotation.
    The result can be depicted by isotools.plots.plot_bar.

    :param groups: A dict {grouname:[sample_name_list]} specifing sample groups. If omitted, the samples are analyzed individually.
    :param weight_by_coverage: If True, each transcript is weighted by the coverage.
    :param min_coverage: Threshold to ignore poorly covered transcripts.
    :param tr_filter: Filter dict, that is passed to self.iter_transcripts().
    :return: Table with numbers of novel alternative splicing events, and suggested parameters for isotools.plots.plot_bar().'''
    weights = dict()
    # if groups is not None:
    #    gi={r:i for i,r in enumerate(runs)}
    #    groups={gn:[gi[r] for r in gr] for gn,gr in groups.items()}
    current = None
    if groups is not None:
        sidx = {sa: i for i, sa in enumerate(self.samples)}  # idx
        groups = {gn: [sidx[sa] for sa in gr] for gn, gr in groups.items()}

    for g, trid, tr in self.iter_transcripts(**tr_filter):
        if g != current:
            current = g
            w = g.coverage.copy() if groups is None else np.array([g.coverage[grp, :].sum(0) for grp in groups.values()])
            w[w < min_coverage] = 0
            if not weight_by_coverage:
                w[w > 0] = 1
        if 'annotation' not in tr or tr['annotation'] is None:
            weights['unknown'] = weights.get('unknown', np.zeros(w.shape[0])) + w[:, trid]
        else:
            for stype in tr['annotation'][1]:
                weights[stype] = weights.get(stype, np.zeros(w.shape[0])) + w[:, trid]
        weights['total'] = weights.get('total', np.zeros(w.shape[0])) + w[:, trid]

    df = pd.DataFrame(weights, index=self.samples if groups is None else groups.keys()).T
    df = df.reindex(df.mean(1).sort_values(ascending=False).index, axis=0)  # sort by row mean
    if weight_by_coverage:
        title = 'Expressed Transcripts'
        ylab = 'fraction of reads'
    else:
        title = 'Different Transcripts'
        ylab = 'fraction of  different transcripts'
        if min_coverage > 1:
            title += f' > {min_coverage} reads'

    return df, {'ylabel': ylab, 'title': title}
    #


def filter_stats(self, tags=None, groups=None, weight_by_coverage=True, min_coverage=2, tr_filter={}):
    '''Summary statistics for filter flags.

    This function counts the number of transcripts correspondign to filter tags.
    The result can be depicted by isotools.plots.plot_bar.

    :param tags: The filter tags to be evaluated. If omitted, all transcript tags are selected.
    :param groups: A dict {grouname:[sample_name_list]} specifing sample groups. If omitted, the samples are analyzed individually.
    :param weight_by_coverage: If True, each transcript is weighted by the number of supporting reads.
    :param min_coverage: Coverage threshold per sample to ignore poorly covered transcripts.
    :param tr_filter: Only transcripts that pass this filter are evaluated. Filter is provided as dict of parameters, passed to self.iter_transcripts().
    :return: Table with numbers of transcripts featuring the filter tag, and suggested parameters for isotools.plots.plot_bar().'''

    weights = dict()
    if tags is None:
        tags = list(self.filter['transcript'])
    assert all(t in self.filter['transcript'] for t in tags), '"Tags" contains invalid tags'
    filterfun = {t: _filter_function(self.filter['transcript'][t])[0] for t in tags}
    if groups is not None:
        sidx = {sa: i for i, sa in enumerate(self.samples)}  # idx
        groups = {gn: [sidx[sa] for sa in gr] for gn, gr in groups.items()}
    current = None
    for g, trid, tr in self.iter_transcripts(**tr_filter):
        if g != current:
            current = g
            w = g.coverage.copy() if groups is None else np.array([g.coverage[grp, :].sum(0) for grp in groups.values()])
            w[w < min_coverage] = 0
            if not weight_by_coverage:
                w[w > 0] = 1
        # relevant_filter=[f for f in tr['filter'] if  consider is None or f in consider]
        relevant_filter = [t for t in tags if filterfun[t](g=g, trid=trid, **tr)]
        for f in relevant_filter:
            weights[f] = weights.get(f, np.zeros(w.shape[0])) + w[:, trid]
        if not relevant_filter:
            weights['PASS'] = weights.get('PASS', np.zeros(w.shape[0])) + w[:, trid]
        weights['total'] = weights.get('total', np.zeros(w.shape[0])) + w[:, trid]

    df = pd.DataFrame(weights, index=self.samples if groups is None else groups.keys()).T

    df = df.reindex(df.mean(1).sort_values(ascending=False).index, axis=0)
    ylab = 'fraction of reads' if weight_by_coverage else 'fraction of different transcripts'
    if weight_by_coverage:
        title = 'Expressed Transcripts'
    else:
        title = 'Different Transcripts'
    if min_coverage > 1:
        title += f' > {min_coverage} reads'
    return df, {'ylabel': ylab, 'title': title}


def transcript_length_hist(self=None, groups=None, add_reference=False, bins=50, x_range=(
        0, 10000), weight_by_coverage=True, min_coverage=2, use_alignment=True, tr_filter={}, ref_filter={}):
    '''Retrieves the transcript length distribution.

    This function counts the number of transcripts within length interval.
    The result can be depicted by isotools.plots.plot_dist.

    :param groups: A dict {grouname:[sample_name_list]} specifing sample groups. If omitted, the samples are analyzed individually.
    :param add_reference: Add the transcript length distribution of the reference annotaiton.
    :param bins: Define the length interval, either by a single number of bins, or by a list of lengths, defining the interval boudaries.
    :param x_range: The range of the intervals. Ignored if "bins" is provided as a list.
    :param weight_by_coverage: If True, each transcript is weighted by the coverage.
    :param min_coverage: Threshold to ignore poorly covered transcripts.
    :param use_alignment: use the transcript length as defined by the alignment (e.g. the sum of all exon lengths).
    :param tr_filter: Filter dict, that is passed to self.iter_transcripts().
    :param ref_filter: Filter dict, that is passed to self.iter_ref_transcripts() (relevant only if add_reference=True).
    :return: Table with numbers of transcripts within the length intervals, and suggested parameters for isotools.plots.plot_distr().'''

    trlen = []
    cov = []
    current = None
    for g, trid, tr in self.iter_transcripts(**tr_filter):
        if g != current:
            current = g
            current_cov = g.coverage
        cov.append(current_cov[:, trid])
        trlen.append(sum(e[1] - e[0] for e in tr['exons']) if use_alignment else tr['source_len'])  # source_len is not set in the current version
    cov = pd.DataFrame(cov, columns=self.samples)
    if groups is not None:
        cov = pd.DataFrame({grn: cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins, int):
        bins = np.linspace(x_range[0] - .5, x_range[1] - .5, bins + 1)
    cov[cov < min_coverage] = 0
    if not weight_by_coverage:
        cov[cov > 0] = 1
    counts = pd.DataFrame({gn: np.histogram(trlen, weights=g_cov, bins=bins)[0] for gn, g_cov in cov.items()})
    if add_reference:
        ref_len = [sum(e[1] - e[0] for e in tr['exons']) for _, _, tr in self.iter_ref_transcripts(**ref_filter)]
        counts['reference'] = np.histogram(ref_len, bins=bins)[0]
    bin_df = pd.DataFrame({'from': bins[:-1], 'to': bins[1:]})
    params = dict(yscale='linear', title='transcript length', xlabel='transcript length [bp]', density=True)
    return pd.concat([bin_df, counts], axis=1).set_index(['from', 'to']), params


def transcript_coverage_hist(self, groups=None, bins=50, x_range=(1, 1001), tr_filter={}):
    '''Retrieves the transcript coverage distribution.

    This function counts the number of transcripts within coverage interval.
    The result can be depicted by isotools.plots.plot_dist.

    :param groups: A dict {grouname:[sample_name_list]} specifing sample groups. If omitted, the samples are analyzed individually.
    :param bins: Define the covarge interval, either by a single number of bins, or by a list of values, defining the interval boudaries.
    :param x_range: The range of the intervals. Ignored if "bins" is provided as a list.
    :param tr_filter: Filter dict, that is passed to self.iter_transcripts().
    :return: Table with numbers of transcripts within the coverage intervals, and suggested parameters for isotools.plots.plot_distr().'''
    # get the transcript coverage in bins for groups
    # return count dataframe and suggested default parameters for plot_distr
    cov = []
    current = None
    for g, trid, _ in self.iter_transcripts(**tr_filter):
        if g != current:
            current = g
            current_cov = g.coverage
        cov.append(current_cov[:, trid])
    cov = pd.DataFrame(cov, columns=self.samples)
    if groups is not None:
        cov = pd.DataFrame({grn: cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins, int):
        bins = np.linspace(x_range[0] - .5, x_range[1] - .5, bins + 1)
    counts = pd.DataFrame({gn: np.histogram(g_cov, bins=bins)[0] for gn, g_cov in cov.items()})
    bin_df = pd.DataFrame({'from': bins[:-1], 'to': bins[1:]})
    params = dict(yscale='log', title='transcript coverage', xlabel='reads per transcript')
    return pd.concat([bin_df, counts], axis=1).set_index(['from', 'to']), params
    # plot histogram
    # cov.mask(cov.lt(x_range[0]) | cov.gt(x_range[1])).plot.hist(ax=ax, alpha=0.5, bins=n_bins)
    # ax=counts.plot.bar()
    # ax.plot(x, counts)


def transcripts_per_gene_hist(self, groups=None, add_reference=False, bins=49, x_range=(1, 50), min_coverage=2, tr_filter={}, ref_filter={}):
    '''Retrieves the histogram of number of transcripts per gene.

    This function counts the genes featuring transcript numbers within specified intervals.
    The result can be depicted by isotools.plots.plot_dist.

    :param groups: A dict {grouname:[sample_name_list]} specifing sample groups. If omitted, the samples are analyzed individually.
    :param add_reference: Add the transcript per gene histogram of the reference annotaiton.
    :param bins: Define the intervals, either by a single number of bins, or by a list of values, defining the interval boudaries.
    :param x_range: The range of the intervals. Ignored if "bins" is provided as a list.
    :param min_coverage: Threshold to ignore poorly covered transcripts.
    :param tr_filter: Filter dict, that is passed to self.iter_transcripts().
    :param ref_filter: Filter dict, that is passed to self.iter_ref_transcripts() (relevant only if add_reference=True).
    :return: Table with numbers of genes featuring transcript numbers within the specified intervals,
        and suggested parameters for isotools.plots.plot_distr().'''
    ntr = []
    current = None
    if groups is None:
        group_names = self.samples
    else:
        group_names = groups.keys()
        sidx = {sa: i for i, sa in enumerate(self.samples)}  # idx
        groups = {gn: [sidx[sa] for sa in gr] for gn, gr in groups.items()}
    n_sa = len(group_names)
    for g, trid, _ in self.iter_transcripts(**tr_filter):
        if g != current:
            current = g
            current_cov = g.coverage if groups is None else np.array([g.coverage[grp, :].sum(0) for grp in groups.values()])
            ntr.append(np.zeros(n_sa))
        ntr[-1] += current_cov[:, trid] >= min_coverage

    ntr = pd.DataFrame((n for n in ntr if n.sum() > 0), columns=group_names)
    if isinstance(bins, int):
        bins = np.linspace(x_range[0] - .5, x_range[1] - .5, bins + 1)
    counts = pd.DataFrame({gn: np.histogram(n, bins=bins)[0] for gn, n in ntr.items()})
    if add_reference:
        if ref_filter:
            logger.warning('reference filter not implemented')
        ref_ntr = [g.n_ref_transcripts for g in self]  # todo: add reference filter
        counts['reference'] = np.histogram(ref_ntr, bins=bins)[0]
    bin_df = pd.DataFrame({'from': bins[:-1], 'to': bins[1:]})
    sub = f'counting transcripts covered by >= {min_coverage} reads'
    if 'query' in tr_filter:
        sub += f', filter query: {tr_filter["query"]}'
    params = dict(yscale='log', title='transcript per gene\n' + sub, xlabel='transcript per gene')
    return pd.concat([bin_df, counts], axis=1).set_index(['from', 'to']), params


def exons_per_transcript_hist(self, groups=None, add_reference=False, bins=34, x_range=(1, 69),
                              weight_by_coverage=True, min_coverage=2, tr_filter={}, ref_filter={}):
    '''Retrieves the histogram of number of exons per transcript.

    This function counts the transcripts featuring exon numbers within specified intervals.
    The result can be depicted by isotools.plots.plot_dist.

    :param groups: A dict {grouname:[sample_name_list]} specifing sample groups. If omitted, the samples are analyzed individually.
    :param add_reference: Add the exons per transcript histogram of the reference annotaiton.
    :param bins: Define the intervals, either by a single number of bins, or by a list of values, defining the interval boudaries.
    :param x_range: The range of the intervals. Ignored if "bins" is provided as a list.
    :param weight_by_coverage: If True, each transcript is weighted by the coverage.
    :param min_coverage: Threshold to ignore poorly covered transcripts.
    :param tr_filter: Filter dict, that is passed to self.iter_transcripts().
    :param ref_filter: Filter dict, that is passed to self.iter_ref_transcripts() (relevant only if add_reference=True).
    :return: Table with numbers of transcripts featuring exon numbers within the specified intervals,
        and suggested parameters for isotools.plots.plot_distr().'''
    n_exons = []
    cov = []
    current = None
    for g, trid, tr in self.iter_transcripts(**tr_filter):
        if g != current:
            current = g
            current_cov = g.coverage
        cov.append(current_cov[:, trid])
        n_exons.append(len(tr['exons']))
    cov = pd.DataFrame(cov, columns=self.samples)
    if groups is not None:
        cov = pd.DataFrame({grn: cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins, int):
        bins = np.linspace(x_range[0] - .5, x_range[1] - .5, bins + 1)
    cov[cov < min_coverage] = 0
    if not weight_by_coverage:
        cov[cov > 0] = 1
    counts = pd.DataFrame({gn: np.histogram(n_exons, weights=g_cov, bins=bins)[0] for gn, g_cov in cov.items()})
    if add_reference:
        ref_n_exons = [len(tr['exons']) for _, _, tr in self.iter_ref_transcripts(**ref_filter)]
        counts['reference'] = np.histogram(ref_n_exons, bins=bins)[0]
    bin_df = pd.DataFrame({'from': bins[:-1], 'to': bins[1:]})
    sub = f'counting transcripts covered by >= {min_coverage} reads'
    if 'query' in tr_filter:
        sub += f', filter query: {tr_filter["query"]}'
    params = dict(yscale='log', title='exons per transcript\n' + sub, xlabel='number of exons per transcript')
    return pd.concat([bin_df, counts], axis=1).set_index(['from', 'to']), params


def downstream_a_hist(self, groups=None, add_reference=False, bins=30, x_range=(0, 1), weight_by_coverage=True, min_coverage=2, tr_filter={}, ref_filter={}):
    '''Retrieves the distribution of downstream adenosine content.

    High downstream adenosine content is indicative for internal priming.

    :param groups: A dict {grouname:[sample_name_list]} specifing sample groups. If omitted, the samples are analyzed individually.
    :param add_reference: Add the distribution of downstream adenosine content of the reference annotaiton.
    :param bins: Define the intervals, either by a single number of bins, or by a list of values, defining the interval boudaries.
    :param x_range: The range of the intervals. Ignored if "bins" is provided as a list. Should not exceed (0,1), e.g. 0 to 100%.
    :param weight_by_coverage: If True, each transcript is weighted by the coverage.
    :param min_coverage: Threshold to ignore poorly covered transcripts.
    :param tr_filter: Filter dict, that is passed to self.iter_transcripts().
    :param ref_filter: Filter dict, that is passed to self.iter_ref_transcripts() (relevant only if add_reference=True).
    :return: Table with downstream adenosine content distribution, and suggested parameters for isotools.plots.plot_distr().'''
    acontent = []
    cov = []
    current = None
    for g, trid, tr in self.iter_transcripts(**tr_filter):
        if g != current:
            current = g
            current_cov = g.coverage
        cov.append(current_cov[:, trid])
        try:
            acontent.append(tr['downstream_A_content'])
        except KeyError:
            acontent.append(-1)
    cov = pd.DataFrame(cov, columns=self.samples)
    if groups is not None:
        cov = pd.DataFrame({grn: cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins, int):
        bins = np.linspace(x_range[0], x_range[1], bins + 1)
    cov[cov <= min_coverage] = 0
    if not weight_by_coverage:
        cov[cov > 0] = 1
    counts = pd.DataFrame({gn: np.histogram(acontent, weights=g_cov, bins=bins)[0] for gn, g_cov in cov.items()})
    if add_reference:
        ref_acontent = [tr['downstream_A_content'] for _, _, tr in self.iter_ref_transcripts(**ref_filter) if 'downstream_A_content' in tr]
        counts['reference'] = np.histogram(ref_acontent, bins=bins)[0]
    bin_df = pd.DataFrame({'from': bins[:-1], 'to': bins[1:]})
    params = dict(title='downstream genomic A content', xlabel='fraction of A downstream the transcript')
    return pd.concat([bin_df, counts], axis=1).set_index(['from', 'to']), params


def direct_repeat_hist(self, groups=None, bins=10, x_range=(0, 10), weight_by_coverage=True, min_coverage=2, tr_filter={}):
    '''Retrieves the distribution direct repeat length at splice junctions.

    Direct repeats are indicative for reverse transcriptase template switching.

    :param groups: A dict {grouname:[sample_name_list]} specifing sample groups. If omitted, the samples are analyzed individually.
    :param bins: Define the intervals, either by a single number of bins, or by a list of values, defining the interval boudaries.
    :param x_range: The range of the intervals. Ignored if "bins" is provided as a list.
    :param weight_by_coverage: If True, each transcript is weighted by the coverage.
    :param min_coverage: Threshold to ignore poorly covered transcripts.
    :param tr_filter: Filter dict, that is passed to self.iter_transcripts().
    :return: Table with direct repeat length distribution, and suggested parameters for isotools.plots.plot_distr().'''
    # find the direct repeat length distribution in FSM transcripts and putative RTTS
    # putative RTTS are identified by introns where both splice sites are novel but within annotated exons
    # TODO: actually no need to check annotation, could simply use filter flags (or the definition from the filter flags, which should be faster)
    rl = {cat: [] for cat in ('known', 'novel canonical', 'novel noncanonical')}
    for g, trid, tr in self.iter_transcripts(**tr_filter):
        if 'annotation' in tr and tr['annotation'][0] == 0:  # e.g. FSM
            rl['known'].extend((drl, g.coverage[:, trid]) for drl in tr['direct_repeat_len'])
        elif g.is_annotated and 'novel_splice_sites' in tr:
            novel_junction = [i // 2 for i in tr['novel_splice_sites'] if i % 2 == 0 and i + 1 in tr['novel_splice_sites']]
            nc = {v[0] for v in tr.get('noncanonical_splicing', [])}
            rl['novel noncanonical'].extend((tr['direct_repeat_len'][sj], g.coverage[:, trid]) for sj in novel_junction if sj in nc)
            rl['novel canonical'].extend((tr['direct_repeat_len'][sj], g.coverage[:, trid]) for sj in novel_junction if sj not in nc)

    rl_cov = {cat: pd.DataFrame((v[1] for v in rl[cat]), columns=self.samples) for cat in rl}
    if groups is not None:
        rl_cov = {cat: pd.DataFrame({grn: rl_cov[cat][grp].sum(1) for grn, grp in groups.items()}) for cat in rl_cov}
    for cov_df in rl_cov.values():
        cov_df[cov_df < min_coverage] = 0
        if not weight_by_coverage:
            cov_df[cov_df > 0] = 1
    if isinstance(bins, int):
        bins = np.linspace(x_range[0] - .5, x_range[1] - .5, bins + 1)
    counts = pd.DataFrame({f'{sa} {cat}': np.histogram([val[0] for val in rl_list], weights=rl_cov[cat][sa], bins=bins)[
                          0] for cat, rl_list in rl.items() for sa in (self.samples if groups is None else groups)})

    bin_df = pd.DataFrame({'from': bins[:-1], 'to': bins[1:]})
    params = dict(title='direct repeat length', xlabel='length of direct repeats at splice junctons', ylabel='# transcripts')

    return pd.concat([bin_df, counts], axis=1).set_index(['from', 'to']), params


def rarefaction(self, groups=None, fractions=20, min_coverage=2, tr_filter={}):
    '''Rarefaction analysis

    Reads are subsampled according to the provided fractions, to estimate saturation of the transcriptome.

    :param groups: A dict {grouname:[sample_name_list]} specifing sample groups. If omitted, the samples are analyzed individually.
    :param fractions: the fractions of reads to be subsampled.
        Either a list of floats between 0 and 1, or a integer number, specifying the number of equally spaced fractions.
    :param min_coverage: Number of reads per transcript required to consider the transcriped discovered.
    :param tr_filter: Filter dict, that is passed to self.iter_transcripts().
    :return: Tuple with:
        1) Data frame containing the number of discovered transcripts, for each subsampling fraction and each sample / sample group.
        2) Dict with total number of reads for each group. '''

    cov = []
    current = None
    if isinstance(fractions, int):
        fractions = np.linspace(1 / fractions, 1, fractions)
    for g, trid, _ in self.iter_transcripts(**tr_filter):
        if g != current:
            current = g
            current_cov = g.coverage
        cov.append(current_cov[:, trid])
    cov = pd.DataFrame(cov, columns=self.samples)
    total = dict(self.sample_table.set_index('name').nonchimeric_reads)
    if groups is not None:
        cov = pd.DataFrame({grn: cov[grp].sum(1) for grn, grp in groups.items()})
        total = {grn: sum(n for sa, n in total.items() if sa in grp) for grn, grp in groups.items()}
    curves = {}
    for sa in cov:
        curves[sa] = [(np.random.binomial(n=cov[sa], p=th) >= min_coverage).sum() for th in fractions]
    return pd.DataFrame(curves, index=fractions), total
