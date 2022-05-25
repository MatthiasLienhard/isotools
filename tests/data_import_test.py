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
    isoseq_csv._add_novel_gene('nix',10,20,'-', {'exons':[10,20]}) #make it a little harder
    id_map=isoseq_csv.add_sample_from_csv('tests/data/example_1_cov.csv', 'tests/data/example_1.gtf', reconstruct_genes=True)
    remapped_genes={gid:gid2 for gid2,id_dict in id_map.items() for gid in id_dict.values()}
    logger.info('remapped %s transcripts', sum(len(d) for d in id_map))
    discrepancy = False
    for g in isoseq.iter_genes(query='EXPRESSED'):
        if (g.is_annotated and g.id in remapped_genes) or (g.id not in isoseq_csv and g.id not in remapped_genes):
            logger.error('gene missing/renamed after csv import: %s' % str(g))
            discrepancy = True
    for g_csv in isoseq_csv.iter_genes(query='EXPRESSED'):
        if not g_csv.is_annotated and g_csv.id in remapped_genes:
            gene_id=remapped_genes[g_csv.id]
        else:
            gene_id=g_csv.id
        g = isoseq[gene_id]
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
    isoseq_csv._add_novel_gene('nix',10,20,'-', {'exons':[10,20]}) #make it a little harder
    id_map=isoseq_csv.add_sample_from_csv('tests/data/example_1_cov.csv', 'tests/data/example_1.gtf', reconstruct_genes=False)
    remapped_genes={v:k for k,v in id_map.items()}
    logger.info('remapped %s genes', len(id_map))
    discrepancy = False
    for g in isoseq.iter_genes(query='EXPRESSED'):
        if (g.is_annotated and g.id in remapped_genes) or (g.id not in isoseq_csv and g.id not in remapped_genes):
            logger.error('gene missing/renamed after csv import: %s' % str(g))
            discrepancy = True
    for g_csv in isoseq_csv.iter_genes(query='EXPRESSED'):
        if not g_csv.is_annotated and g_csv.id in remapped_genes:
            gene_id=remapped_genes[g_csv.id]
        else:
            gene_id=g_csv.id
        g = isoseq[gene_id]
        if len(g.transcripts) != len(g_csv.transcripts):
            logger.error('number of transcripts for %s changed after csv import: %s != %s', g.id, len(g.transcripts), len(g_csv.transcripts))
            discrepancy = True
    assert not discrepancy, 'discrepancy found after csv import'
