import pytest
from isotools.transcriptome import Transcriptome
from isotools._utils import _find_splice_sites, _get_overlap, _get_exonic_region

# @pytest.mark.dependency(depends=['test_import_bam'])


def test_import_find_splice_site():
    isoseq = Transcriptome.load('tests/data/example_1_isotools.pkl')
    for g, _, tr in isoseq.iter_transcripts(query='not NOVEL_GENE'):
        c1 = g.ref_segment_graph.find_splice_sites(tr['exons'])
        c2 = _find_splice_sites(tr['exons'], g.ref_transcripts)
        assert all(c1 == c2), 'isotools._transcriptome_io._find_splice_sites and Segment_Graph.find_splice_sites yield different results'

@pytest.mark.dependency()
def test_exon_regions():
    isoseq = Transcriptome.load('tests/data/example_1_isotools.pkl')
    for g in isoseq.iter_genes(query='not NOVEL_GENE'):
        c1 = g.ref_segment_graph.get_exonic_region()
        c2 = _get_exonic_region(g.ref_transcripts)
        assert len(c1) == len(c2), 'isotools._transcriptome_io._get_exonic_region and Segment_Graph.get_exonic_region yield different length'
        assert all(reg1[0] == reg2[0] and reg1[1] == reg2[1] for reg1, reg2 in zip(c1, c2)),\
            'isotools._transcriptome_io._get_exonic_region and Segment_Graph.get_exonic_region yield different regions'
    assert True


@pytest.mark.dependency(depends=['test_exon_regions'])
def test_import_exonic_overlap():
    isoseq = Transcriptome.load('tests/data/example_1_isotools.pkl')
    for g, _, tr in isoseq.iter_transcripts(query='not NOVEL_GENE'):
        c1 = g.ref_segment_graph.get_overlap(tr['exons'])[0]
        c2 = _get_overlap(tr['exons'], g.ref_transcripts)
        assert c1 == c2, 'isotools._transcriptome_io._get_overlap and Segment_Graph.get_overlap yield different results'
    assert True