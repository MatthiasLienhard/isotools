import isotools.transcriptome
import unittest

class TestAltsplice(unittest.TestCase):
    def setUp(self):
        tr1={'exons':[[21,30], [41,50], [61,65]]}
        tr2={'exons':[[21,30], [41,50], [71,75]]}
        self.rg1=isotools.transcriptome.Gene(21,85,{'transcripts':[tr1,tr2],'chr':'chr1','strand':'+', 'Name':'RefGene1','ID':'R1'})
        
    def test_id(self):
        exons=[[15,30], [41,50], [61,66]]
        cat,alt_splice=self.rg1.splice_graph.get_alternative_splicing(exons, '+')
        self.assertEqual(cat,0, 'test should return category 0')
        self.assertEqual(alt_splice,{'FSM':0}, 'test find FSM with transcript 0')
'''
    def test_exon_skipping(self):
        exons=[[15,30],  [61,66]]
        alt_splice=isotools.transcriptome.get_splice_type(self.ref,exons, False)
        self.assertEqual(len(alt_splice),1, 'test should find exactly one alternative splicing event')
        self.assertEqual(alt_splice, {'exon_skipping': ['30~61']})

    def test_gapped_exon(self):
        exons=[[21,23],[26,30],[41,50], [61,70]]
        alt_splice=isotools.transcriptome.get_splice_type(self.ref, exons, False)
        self.assertEqual(len(alt_splice),1, 'test should find exactly one alternative splicing event')
        assert 'gapped_exon' in alt_splice,'gapped exon not found'
        self.assertEqual(alt_splice['gapped_exon'][0],'21-23~26-30','position of gapped exon is wrong')

    def test_truncations(self):
        exons=[ [45,50], [61,70]]
        alt_splice=isotools.transcriptome.get_splice_type(self.ref, exons, False)
        self.assertEqual(alt_splice,{'truncated5': ['45-50']}, '5 truncation at lower end not found')
        alt_splice=isotools.transcriptome.get_splice_type(self.ref, exons, True)
        self.assertEqual(alt_splice,{'truncated3': ['45-50']}, '3 truncation at lower end not found')
        exons=[[21,30], [41,50]]
        alt_splice=isotools.transcriptome.get_splice_type(self.ref, exons, False)
        self.assertEqual(alt_splice,{'truncated3': ['41-50']}, '3 truncation at higer end not found')

    def test_novel_exon(self):
        exons=[[21,30], [41,50],[53,56], [61,70]]
        alt_splice=isotools.transcriptome.get_splice_type(self.ref, exons, False)
        self.assertEqual(alt_splice,{'novel_exon': ['53-56']}, 'novel exon not found')

    def test_novel_tss(self): #'alternative_promoter', 'alternative_polyA'
        exons=[ [35,50], [61,70]]
        alt_splice=isotools.transcriptome.get_splice_type(self.ref, exons, False)
        self.assertEqual(alt_splice,{'alternative_promoter': ['35-50']}, 'alternative promoter not found')
        exons=[[21,30], [41,55]]
        alt_splice=isotools.transcriptome.get_splice_type(self.ref, exons, False)
        self.assertEqual(alt_splice,{'alternative_polyA': ['41-55']}, 'alternative polyA site not found')
        

    def test_retained_intron(self):
        exons=[[21,50], [61,70]]
        alt_splice=isotools.transcriptome.get_splice_type(self.ref, exons, False)
        #fails (also reports 'alternative_donor': [(21, -20)], 'alternative_acceptor': [(50, 20)],)
        self.assertEqual(alt_splice,{'retained_intron': ['21-50']}, 'retained intron not found')
        exons=[[24,55], [61,70]]
        alt_splice=isotools.transcriptome.get_splice_type(self.ref, exons, False)
        #this should report the alternative acceptor
        self.assertEqual(alt_splice,{ 'alternative_acceptor': [(55, 5)], 'retained_intron': ['24-55']},'retained intron with alternative acceptor not found')
        


    def test_alternative_acceptor(self):
        exons=[[15,35], [41,50], [61,66]]
        alt_splice=isotools.transcriptome.get_splice_type(self.ref, exons, False)
        self.assertEqual(alt_splice,{'alternative_acceptor': [(35, 5)]}, 'alternative acceptor not found')
    '''
if __name__ == "__main__":
    unittest.main()

import sys
from importlib import reload

sys.path.append('/project/42/pacbio/bin/isotools/')
import  isotools.splice_graph
reload(isotools.splice_graph)
import logging

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

tr1=[[21,30], [61,65], [71,80]]
tr2=[[21,30], [41,50], [71,75]]
reload(isotools.splice_graph)
sg=isotools.splice_graph.SpliceGraph([tr1,tr2])
exons=[[15,30], [41,50], [71,76]]
cat,alt_splice=sg.get_alternative_splicing(exons, '+')
print(f'FSM: {cat} {alt_splice}')

exons=[[15,30],  [71,86]]
cat,alt_splice=sg.get_alternative_splicing(exons, '+')
print(f'exon skipping: {cat} {alt_splice}')

exons=[[15,32],  [71,86]]
cat,alt_splice=sg.get_alternative_splicing(exons, '+')
print(f'exon skipping with alt donor: {cat} {alt_splice}')

exons=[[15,50]]
cat,alt_splice=sg.get_alternative_splicing(exons, '+')
print(f'mono-exon with retention: {cat} {alt_splice}')

exons=[[21,23],[26,30],[41,50], [71,80]]
cat,alt_splice=sg.get_alternative_splicing(exons, '+')
print(f'gapped_exon: {cat} {alt_splice}')

exons=[ [45,50], [71,80]]
cat,alt_splice=sg.get_alternative_splicing(exons, '+')
print(f'fraction5\': {cat} {alt_splice}')

exons=[[21,30], [41,50]]
cat,alt_splice=sg.get_alternative_splicing(exons, '+')
print(f'fraction3\': {cat} {alt_splice}')

exons=[[21,30], [41,51]]
cat,alt_splice=sg.get_alternative_splicing(exons, '+')
print(f'fraction3 als alt PAS\': {cat} {alt_splice}')


exons=[[21,30], [41,50], [53,56], [71,80]]
cat,alt_splice=sg.get_alternative_splicing(exons, '+')
print(f'novel exon: {cat} {alt_splice}')


exons=[ [35,50], [71,80]]
cat,alt_splice=sg.get_alternative_splicing(exons, '+')
print(f'novel TSS: {cat} {alt_splice}')
exons=[[21,30], [41,55]]
cat,alt_splice=sg.get_alternative_splicing(exons, '+')
print(f'Novel PAS: {cat} {alt_splice}')

exons=[[21,50], [71,80]]
cat,alt_splice=sg.get_alternative_splicing(exons, '+')
print(f'retained intron: {cat} {alt_splice}')
exons=[[24,55], [71,80]]
cat,alt_splice=sg.get_alternative_splicing(exons, '+')
print(f'retained intron with alt donor: {cat} {alt_splice}')

exons=[[15,35], [41,50], [61,66]]
cat,alt_splice=sg.get_alternative_splicing(exons, '+')
print(f'alternative donor: {cat} {alt_splice}')

