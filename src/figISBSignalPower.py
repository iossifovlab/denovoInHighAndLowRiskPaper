#!/usr/bin/env python

from pylab import *
from scipy.stats import poisson,binom_test

I = 1000

BRS = [37529/1865.,3002/1865.]
EFS = [0.05, 0.1, 0.2, 0.3, 0.4]
NS = [1865]
AS = [0.05]

ALT = ['-', '.-', '*-']

R = {}
for bri,br in enumerate(BRS):
    for efi,ef in enumerate(EFS):
        for ni,n in enumerate(NS):
            Ps = poisson.rvs(n*(br+ef),size=I)
            Ss = poisson.rvs(n*br,size=I)
            ps = array([binom_test([p,s]) for p,s in zip(Ps,Ss)])
            Q = [(ps<a).sum() for a in AS]
            R[br,ef,n] = Q

clf()

contribution = array(EFS)*100

for br,lbl in zip(BRS,['all genes','autism LGD target genes']):
    pwrs = [R[br,e,NS[0]][0]/float(I) for e in EFS]
    plot(contribution,pwrs,'*-',label=("%s" % (lbl)))
xticks(contribution,map(lambda c: "%.0f%%" % (c),contribution))
xlabel('contribution')
ylabel('power in %d quads' % (NS[0]))
legend()
gcf().set_size_inches(6,3)
tight_layout()
if len(sys.argv) > 1:
    gcf().savefig(sys.argv[1])
else:
    show()
