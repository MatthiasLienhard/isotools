from isotools import Transcriptome
from isotools.domains import add_domains_to_table


def test_anno_domains():
    isoseq = Transcriptome.load('tests/data/example_1_isotools.pkl')
    isoseq.add_annotation_domains('tests/data/example_anno_domains.csv', category='domains',  progress_bar=False)
    g = isoseq['FN1']
    ref_tr = next(tr for tr in g.ref_transcripts if tr['transcript_name'] == 'FN1-207')
    dom = ref_tr['domain']['annotation']
    assert len(dom) == 25, 'expected 25 annotation domains in FN1-207, but found {len(dom)}'

    diffexpr = isoseq.altsplice_test(groups=isoseq.groups()).sort_values('pvalue')
    diffexpr = add_domains_to_table(diffexpr, isoseq, 'annotation', insert_after='nmdB')
