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
        [[21,30], [37,50], [61,70],[81,90]],    #alternative 3'
        [[21,30], [81,90]]]
weights=np.array([[10,10,10,10,10,10],[5,20,0,15,11,9]])
sg=isotools.splice_graph.SpliceGraph(exons, weights=weights)
list(sg.get_splice_coverage())
#list(sg.get_splice_coverage_old())


[a+b for a in "12" for b in "AB"]
#['1A', '1B', '2A', '2B']
x1=np.array([176.0, 121.0, 87.0, 119.0, 77.0, 188.0, 82.0, 113.0, 210.0, 225.0])
x2=np.array([23.0, 51.0, 240.0, 88.0, 27.0, 45.0, 110.0, 115.0, 11.0])
n1=np.array([200.0, 143.0, 104.0, 139.0, 88.0, 210.0, 96.0, 126.0, 244.0, 249.0])
n2=np.array([83.0, 142.0, 265.0, 218.0, 142.0, 136.0, 165.0, 164.0, 38.0])
res2=(0.4761144448906714, 0.19394337282466823) 	
res1=(0.873307668631287, 0.0014929882277201502)
res0=(0.6652330116312957, 0.230021439336415)
​
import isotools.stats
from scipy.optimize import minimize

reload(isotools.stats)
isotools.stats.betabinom_lr_test([x1,x2],[n1,n2])


p,params=isotools.stats.betabinom_lr_test([x1,x2],[n1,n2])
params_alt=[((-m*(m**2-m+v))/v,((m-1)*(m**2-m+v))/v)  for m,v in params]
(x2/n2).var()-params[1][1] #biological variance
(np.random.binomial(n2.astype(int),params[1][0],(1000,len(n2)))/n2).var() #technical variance



np.random.binomial(n2,params[1][0])

prob=x1/n1
m=prob.mean()
d=max(prob.var(),1e-6) 
e=(m**2-m+d) #helper          

mle = minimize(isotools.stats.loglike_betabinom, x0=[-m*e/d,((m-1)*e)/d],bounds=((1e-6,None),(1e-6,None)),  args=(x1,n1),options={'maxiter': 250})
mle = minimize(isotools.stats.loglike_betabinom, x0=[-m*e/d,((m-1)*e)/d],bounds=((1e-6,None),(1e-6,None)),  args=(x1,n1),options={'maxiter': 250}, jac=isotools.stats.loglike_betabinom_jac)

isotools.stats.loglike_betabinom_alt(params[0],x1,n1)
isotools.stats.loglike_betabinom(params_alt,x1,n1)
mle = minimize(isotools.stats.loglike_betabinom_alt, x0=[minit,dinit],bounds=((1e-6,1-1e-6),(1e-6,None)),  args=(x1,n1),options={'maxiter': 250})

mle = minimize(isotools.stats.loglike_betabinom2, x0=[-d/(m*e),d/((m-1)*e)],bounds=((1e-6,None),(1e-6,None)),  args=(x1,n1),options={'maxiter': 250})

isotools.stats.loglike_betabinom_alt(res1,x1,n1)

ni=n2.values
xi=x2.values

#https://en.wikipedia.org/wiki/Beta-binomial_distribution point estimates, method of moments
m1=(xi/ni)
m2=(xi**2/ni)
a=(ni*m1-m2)/(ni*(m2/m1-m1-1)+m1)
b=(ni-m1)*(ni-m2/m1)/((ni*m2/m1-m1-1)+m1)

mu=a/(a+b)
th=1/(a+b+1)






from intervaltree import IntervalTree, Interval

class Gene(Interval):
    def __new__(cls,begin, end, data, transcriptome):
        return super().__new__(cls,begin, end, data) #required as Interval (and Gene) are immutable
    def __init__(self,begin, end, data, transcriptome):
        self._transcriptome=transcriptome

Gene(1,2,[],5)


from pysam import TabixFile, AlignmentFile, FastaFile
fn='~/projects/molgen/project42/pacbio/references/gencode/gencode.v29.annotation.gff3.gz'
gff=TabixFile(fn)