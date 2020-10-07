import isotools.transcriptome
import unittest
import logging
#isotools.transcriptome.log.setLevel(logging.DEBUG)
from intervaltree import IntervalTree
FORMAT = '%(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger('unit test')
logger.setLevel(logging.INFO)



class TestAnnotation(unittest.TestCase):
    def setUp(self):
        rg1=isotools.transcriptome.Gene(21,65,{'transcripts':{'reftr1.1':{'exons':[[21,30], [41,50], [61,65]]}},'chr':'chr1','strand':'+', 'Name':'RefGene1','ID':'R1'})
        rg2=isotools.transcriptome.Gene(21,65,{'transcripts':{'reftr2.1':{'exons':[[21,30], [41,50], [61,65]]}},'chr':'chr2','strand':'-', 'Name':'RefGene2','ID':'R2'})
        pg1=isotools.transcriptome.Gene(21,65,{'transcripts':{
            'pb_tr1.1':{'exons':[[11,35], [41,50], [61,65]],'fusion':[{'chr':'chr3', 'exons':[[21,30], [41,50], [61,65]],'strand':'-'}]} ,
            'pb_tr1.2':{'exons':[[21,30], [61,65]]},
            'pb_tr1.3':{'exons':[[21,30], [41,50], [61,65]], 'fusion':[{'chr':'chr2', 'exons':[[21,30], [41,50], [61,65]],'strand':'-'}]}},
            'chr':'chr1','strand':'+', 'Name':'PB_Gene1','ID':'PB1'})
        pg2=isotools.transcriptome.Gene(21,65,{'transcripts':{'pb_tr2.1':{'exons':[[25,30], [41,50], [61,70]]}},'chr':'chr2','strand':'-', 'Name':'PB_Gene2','ID':'PB2'})
        pg3=isotools.transcriptome.Gene(21,65,{'transcripts':{'pb_tr2.1':{'exons':[[25,30], [41,50], [61,70]]}},'chr':'chr3','strand':'-', 'Name':'PB_Gene3','ID':'PB3'})#novel gene
        reference=isotools.transcriptome.Transcriptome(data={'chr1':IntervalTree([rg1]),'chr2':IntervalTree([rg2])},infos={'file_name':'NA'})
        self.isoseq=isotools.transcriptome.Transcriptome(data={'chr1':IntervalTree([pg1]),'chr2':IntervalTree([pg3]),'chr3':IntervalTree([pg3])},infos={'file_name':'NA'})
        self.isoseq.annotate(reference)
    
    def test_len(self):
        self.assertEqual(len(self.isoseq)==3, 'len is %i, should be 3' % len(self.isoseq))

    def test_nameing(self):
        self.assertIn(  'R1' , self.isoseq._idx, 'annotation gene naming issue with R1')

    def test_transcript_nameing(self):
        self.assertIn(  'pb_tr1.3' , self.isoseq['R1'].transcripts, 'did not find transcript pb_tr1.3 which should be in gene R1')

    def test_fusion(self): 
        self.assertIn( 'fusion' , self.isoseq['R1'].transcripts['pb_tr1.3'], 'pb_tr1.3 does not have fusions anymore')
        self.assertIn( 'annotation' , self.isoseq['R1'].transcripts['pb_tr1.3']['fusion'][0], 'fusion not annotated')
        f_refgene=self.isoseq['R1'].transcripts['pb_tr1.3']['fusion'][0]['annotation']['ref_gene_name']
        self.assertEqual( f_refgene ,'RefGene2', 'fusion annotation is wrong ("%s" != "RefGene2"' % f_refgene)
        logger.info('fusion transcript: %s' % str(self.isoseq["R1"].transcripts["pb_tr1.1"]))


if __name__ == "__main__":
    unittest.main()




