import pytest
from isotools.transcriptome import Transcriptome
from isotools._utils import splice_identical


@pytest.mark.dependency()
def test_import_gff():
    isoseq = Transcriptome.from_reference('tests/data/example.gff.gz')
    assert len(isoseq)==39, 'we expect 39 genes'
    isoseq.save_reference('tests/data/example_ref_isotools.pkl')
    assert True


@pytest.mark.dependency(depends=['test_import_gff'])
def test_import_bam():
    isoseq = Transcriptome.from_reference('tests/data/example_ref_isotools.pkl')
    assert isoseq.n_transcripts==0, 'there should not be any transcripts'
    for sa in ('CTL', 'VPA'):
        isoseq.add_sample_from_bam(f'tests/data/example_1_{sa}.bam', sample_name=sa, group=sa, platform='SequelII')
    assert isoseq.n_transcripts==185, 'we expect 185 transcripts'
    isoseq.add_qc_metrics('tests/data/example.fa')
    isoseq.save('tests/data/example_1_isotools.pkl')


@pytest.mark.dependency(depends=['test_import_bam'])
def test_fsm():
    isoseq = Transcriptome.load('tests/data/example_1_isotools.pkl')    
    for g, _, tr in isoseq.iter_transcripts(query='FSM'):
        assert tr['annotation'][0] == 0
        for ref_id in tr['annotation'][1]['FSM']:
            assert splice_identical(tr['exons'], g.ref_transcripts[ref_id]['exons'])
