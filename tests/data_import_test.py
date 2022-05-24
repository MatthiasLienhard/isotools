import pytest
from isotools.transcriptome import Transcriptome
from isotools._utils import splice_identical
import logging
logger = logging.getLogger('isotools')


@pytest.mark.dependency()
def test_import_gff():
    isoseq = Transcriptome.from_reference('tests/data/example.gff.gz')
    assert len(isoseq) == 65, 'we expect 65 genes'
    isoseq.save_reference('tests/data/example_ref_isotools.pkl')
    assert True


@pytest.mark.dependency(depends=['test_import_gff'])
def test_import_bam():
    isoseq = Transcriptome.from_reference('tests/data/example_ref_isotools.pkl')
    assert isoseq.n_transcripts == 0, 'there should not be any transcripts'
    for sa in ('CTL', 'VPA'):
        isoseq.add_sample_from_bam(f'tests/data/example_1_{sa}.bam', sample_name=sa, group=sa, platform='SequelII')
    # assert isoseq.n_transcripts == 185, 'we expect 185 transcripts'
    isoseq.add_qc_metrics('tests/data/example.fa')
    isoseq.save('tests/data/example_1_isotools.pkl')


@pytest.mark.dependency(depends=['test_import_bam'])
def test_fsm():
    isoseq = Transcriptome.load('tests/data/example_1_isotools.pkl')
    count = 0
    for g, _, tr in isoseq.iter_transcripts(query='FSM'):
        assert tr['annotation'][0] == 0
        count += 1
        for ref_id in tr['annotation'][1]['FSM']:
            assert splice_identical(tr['exons'], g.ref_transcripts[ref_id]['exons'])
    assert count == 22, 'expected 22 FSM transcripts'


@pytest.mark.dependency(depends=['test_import_bam'])
def test_import_csv_reconstruct():  # reconstruct gene structure from scratch
    isoseq = Transcriptome.load('tests/data/example_1_isotools.pkl')
    cov_tab = isoseq.transcript_table(coverage=True)
    cov_tab.to_csv('tests/data/example_1_cov.csv')
    isoseq.write_gtf('tests/data/example_1.gtf')
    isoseq_csv = Transcriptome.from_reference('tests/data/example_ref_isotools.pkl')
    isoseq_csv.add_sample_from_csv('tests/data/example_1_cov.csv', 'tests/data/example_1.gtf')
    discrepancy = False
    for g in isoseq:
        if g.is_expressed and g.ref_transcripts:
            if g.id not in isoseq_csv:
                logger.error('gene missing after csv import: %s' % str(g))
                discrepancy = True
            g_csv = isoseq_csv[g.id]
            if len(g.transcripts) != len(g_csv.transcripts):
                logger.error('number of transcripts for %s changed after csv import: %s != %s', g.id, len(g.transcripts), len(g_csv.transcripts))
                discrepancy = True
    assert not discrepancy, 'discrepancy found after csv import'


@pytest.mark.dependency(depends=['test_import_bam'])
def test_import_csv():  # use gene structure from gtf
    isoseq = Transcriptome.load('tests/data/example_1_isotools.pkl')
    cov_tab = isoseq.transcript_table(coverage=True)
    cov_tab.to_csv('tests/data/example_1_cov.csv')
    isoseq.write_gtf('tests/data/example_1.gtf')
    isoseq_csv = Transcriptome.from_reference('tests/data/example_ref_isotools.pkl')
    isoseq_csv.add_sample_from_csv('tests/data/example_1_cov.csv', 'tests/data/example_1.gtf', reconstruct_genes=False)
    discrepancy = False
    for g in isoseq:
        if g.is_expressed and g.ref_transcripts:
            if g.id not in isoseq_csv:
                logger.error('gene missing after csv import: %s' % str(g))
                discrepancy = True
            g_csv = isoseq_csv[g.id]
            if len(g.transcripts) != len(g_csv.transcripts):
                logger.error('number of transcripts for %s changed after csv import: %s != %s', g.id, len(g.transcripts), len(g_csv.transcripts))
                discrepancy = True
    assert not discrepancy, 'discrepancy found after csv import'
