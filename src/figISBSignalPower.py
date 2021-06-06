#!/usr/bin/env python

from pylab import *
from scipy.stats import poisson,binom_test

outFigName = None if len(sys.argv) < 2 else sys.argv[1]
seedV = None if len(sys.argv) < 3 else int(sys.argv[2])
IFN = "intronic_result_table.txt" if len(sys.argv) < 4 else int(sys.argv[3])

if seedV: seed(seedV)

RS = genfromtxt(IFN, delimiter='\t',dtype=None,names=True, case_sensitive=True)
RSD = {}
for R in RS:
    k = (R['set'],R['eventType'],R['intronType'])
    RSD[k] = R

allR = RSD["all genes",'sub','inter-coding_intronic']
autR = RSD["autism LGD",'sub','inter-coding_intronic']

nIntronSubsUnaffectedAll = allR['sRealNu'] 
nIntronSubsUnaffectedCandidateGenes = autR['sRealNu'] 
nUnaffected = allR['sRealU'] # all reacords have the same number of affected and unaffected children 
nAffected  = allR['sRealA'] 

# background rates
BRS = [nIntronSubsUnaffectedAll/float(nUnaffected),
       nIntronSubsUnaffectedCandidateGenes/float(nUnaffected)]

# effects sizes
EFS = [0.05, 0.1, 0.2, 0.3, 0.4]

# significances
AS = [0.05]

# sample sizes
NS = [nAffected]


I = 1000

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
if outFigName:
    gcf().savefig(outFigName)
else:
    show()
