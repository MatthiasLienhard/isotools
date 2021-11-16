import pytest
from isotools import Gene


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
