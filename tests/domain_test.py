from isotools import Transcriptome


def test_anno_domains():
    isoseq = Transcriptome.load('tests/data/example_1_isotools.pkl')
    isoseq.add_annotation_domains('tests/data/example_anno_domains.csv', category='domains',  progress_bar=False)
    g = isoseq['FN1']
    ref_tr = next(tr for tr in g.ref_transcripts if tr['transcript_name'] == 'FN1-207')
    dom = ref_tr['domain']['annotation']
    assert len(dom) == 25, 'expected 25 annotation domains in FN1-207, but found {len(dom)}'
