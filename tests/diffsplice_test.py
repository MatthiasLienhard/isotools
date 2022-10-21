from isotools import Transcriptome


def test_diffsplice():
    isoseq = Transcriptome.load('tests/data/example_1_isotools.pkl')
    res = isoseq.altsplice_test(groups=isoseq.groups()).sort_values('pvalue')
    assert sum(res.padj < .1) == 1, 'expected exactly one significant event'
    best = res.iloc[0]
    assert best.gene == 'SLC39A14', 'Best differential splicing should be in SLC39A14'
    assert (best.start, best.end) == (408496, 414779), 'Genomic coordinates do not match expectations.'
    assert best.splice_type == 'ME', 'Splice type does not match expectations.'
    assert best.padj < .001, 'Event should be more significant.'
