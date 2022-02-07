import pytest
from isotools import Gene
import numpy as np


@pytest.fixture(scope="session")
def example_gene():
    ref = [[(12, 20), (30, 40), (50, 60), (70, 81)],
           [(11, 20), (35, 40),           (75, 79)],
           [(10, 20), (30, 40), (50, 60), (70, 80)]]
    novel = {'FSM':         [(10, 20), (30, 40), (50, 60), (70, 80)],
             "5' fragment": [(33, 40), (50, 60), (70, 80)],
             "3' fragment": [(10, 20), (30, 40), (50, 55)],
             "mono-exon": [(22, 35)],
             "exon skipping":  [(10, 20), (50, 60), (70, 80)],
             "intron retention":  [(10, 40), (50, 60), (70, 80)],
             "novel combination":  [(10, 20), (30, 40), (75, 80)],
             "novel junction":   [(10, 20), (30, 40), (50, 60), (75, 80)],
             "novel exonic TSS":  [(26, 40), (50, 60), (70, 80)],
             "novel exonic PAS":  [(10, 20), (30, 40), (50, 66)],
             "novel 5' splice site": [(10, 24), (30, 40), (50, 60), (70, 80)],
             "novel 3' splice site": [(10, 20), (26, 40), (50, 60), (70, 80)],
             "novel exon":        [(10, 20), (30, 40), (43, 47), (50, 60), (70, 80)],
             "novel intronic TSS": [(43, 47), (50, 60), (70, 80)],
             "novel intronic PAS": [(10, 20), (30, 40), (82, 90)]}
    ref_tr = [{'exons': e, 'id': f'reference {i+1}'} for i, e in enumerate(ref)]
    transcripts = [{'exons': e, 'transcript_name': n} for n, e in novel.items()]
    g = Gene(10, 81, {'chr': 'chr1', 'strand': '+', 'ID': 'example',
                      'reference': {'transcripts': ref_tr}, 'transcripts': transcripts}, None)
    return g


@pytest.fixture(scope="session")
def example_gene_coor():
    ref = [[(0, 10), (20, 30), (40, 50), (60, 70), (80, 90), (100, 110), (120, 130)]]
    novel = {'priA_priB': [(0, 10), (20, 30), (40, 50), (60, 70), (100, 110), (120, 130)],
             'priA_altB': [(0, 10), (20, 30), (40, 50), (60, 70), (80, 90), (100, 110), (120, 130)],
             'altA_priB': [(0, 10), (20, 50), (60, 70), (100, 110), (120, 130)],
             'altA_altB': [(0, 10), (20, 50), (60, 70), (80, 90), (100, 110), (120, 130)]}

    ref_tr = [{'exons': e, 'id': f'reference {i+1}'} for i, e in enumerate(ref)]
    transcripts = [{'exons': e, 'transcript_name': n} for n, e in novel.items()]
    coverage = np.array([300, 60, 100, 350])
    g = Gene(1, 100, {'chr': 'chr1', 'strand': '+', 'ID': 'example_02',
                      'reference': {'transcripts': ref_tr}, 'transcripts': transcripts,
                      'coverage': coverage}, None)
    return g


@pytest.fixture(scope="session")
def example_gene_uncoor(example_gene_coor):
    g = example_gene_coor
    g.data["ID"] = "example_03"
    g.data['coverage'] = np.array([310, 380, 310, 350])
    return g


@pytest.fixture(scope="session")
def res_coor():
    return [(1.541075314331383e-71, 17.5, 300, 60, 100, 350, 'example_02', 'example_02', 'IR', 'ES', 20, 50, 60, 110)]


@pytest.fixture(scope="session")
def res_uncoor():
    return [(0.47764686057027506, 0.9210526315789473, 310, 380, 310, 350, 'example_03', 'example_03', 'IR', 'ES', 20, 50, 60, 110)]
