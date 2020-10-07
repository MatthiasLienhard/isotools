import isotools.splice_graph
import unittest
from random import sample, randint

#from importlib import reload
#reload(isotools.splice_graph)

class TestSpliceGraph(unittest.TestCase):
    def test_simple(self):
        exons=[[[21,30], [41,50], [61,70]],[[25,30], [41,70], [81,90]],[[25,30], [41,70], [81,90]]]
        self.generic_test(exons)

    def test_random(self):
        n=25
        nsj=400
        sj=sorted(sample(range(1000),nsj))
        exons=[[]]
        for _ in range(n):
            if len(exons[-1])>0:
                exons.append([])
            for j in range(nsj//2):
                r=randint(0,10)
                if r > 8:
                    exons[-1].append([sj[j*2], sj[j*2+1]])
                elif r > 5 and len(exons[-1])>0:
                    exons[-1][-1][1]=sj[j*2+1]
        self.generic_test(exons)


    def generic_test(self, exons):        
        sg=isotools.splice_graph.SpliceGraph(exons)
        for i,tr in enumerate(exons):
            self.assertEqual(tr, sg.restore(i))
            self.assertEqual(tr, sg.restore_reverse(i))

if __name__ == "__main__":
    unittest.main()

