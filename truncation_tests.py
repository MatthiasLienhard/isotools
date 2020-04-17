import isotools.transcriptome
import unittest

class TestSpliceGraph(unittest.TestCase):
    short_tr=[[21,30], [41,50], [61,65]]

    def test_positive(self):
        positive=[[[21,30], [41,50], [61,70]],
                [[21,30], [41,50], [61,70], [81,90]],
                [[25,30], [41,50], [61,70], [81,90]],
                [[19,30], [41,50], [61,70], [81,90]]]
        for long_tr in positive:
            trunc,_,_=isotools.transcriptome.is_truncation(long_tr, self.short_tr)
            self.assertEqual(trunc, True)

    def test_negative(self):
        negative=[[[21,30],[34,36], [41,50], [61,70]],
                [[21,31], [41,50], [61,65]],
                [[21,30], [40,50], [61,65]],
                [[21,30], [41,51], [61,65]]]
        for long_tr in negative:
            trunc,_,_=isotools.transcriptome.is_truncation(long_tr, self.short_tr)
            self.assertEqual(trunc, False)
        

if __name__ == "__main__":
    unittest.main()

