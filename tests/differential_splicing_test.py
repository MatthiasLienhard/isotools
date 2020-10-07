import isotools.splice_graph
import unittest
import numpy as np
from random import sample, randint

class TestSpliceGraph(unittest.TestCase):
    def test_simple(self):
        exons= [[[21,30], [41,50], [61,70],[81,90]], #all exons
                [[21,30], [41,70], [81,90]],    #intron retention
                [[21,30], [61,70], [81,90]],    #exon skipping
                [[21,35], [41,50], [61,70],[81,90]],    #alternative 5'
                [[21,30], [37,50], [61,70],[81,90]]]    #alternative 3'
        weights=np.array([[10,10,10,10,10],[5,20,0,15,10]])
        sg=isotools.splice_graph.SpliceGraph(exons, weights=weights)
        for i,tr in enumerate(exons):
            self.assertEqual(tr, sg.restore(i))
            self.assertEqual(tr, sg.restore_reverse(i))

if __name__ == "__main__":
    unittest.main()

import logging
from importlib import reload
reload(isotools.splice_graph)

isotools.splice_graph.log.setLevel(logging.DEBUG)
exons= [[[21,30], [41,50], [61,70],[81,90]], #all exons
                [[21,30], [41,70], [81,90]],    #intron retention
                [[21,30], [61,70], [81,90]],    #exon skipping
                [[21,35], [41,50], [61,70],[81,90]],    #alternative 5'
                [[21,30], [37,50], [61,70],[81,90]]]    #alternative 3'
weights=np.array([[10,10,10,10,10],[5,20,0,15,10]])
sg=isotools.splice_graph.SpliceGraph(exons, weights=weights)
list(sg.get_splice_coverage())
list(sg.get_splice_coverage_old())

