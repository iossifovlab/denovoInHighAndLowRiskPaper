#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')

from pylab import *
from diData import persons, loadEVS
from result_tables import compareN, empStats

outFigName = None if len(sys.argv) < 2 else sys.argv[1]
seedV = None if len(sys.argv) < 3 else int(sys.argv[2])

if seedV: seed(seedV)

prbs = {pd.pId for pd in persons.values() if pd.coll == 'SSC' and pd.role == 'prb'}
sibs = {pd.pId for pd in persons.values() if pd.coll == 'SSC' and pd.role == 'sib'}

CNVs = loadEVS(['De novo CNV in SSC and AGRE'])

def countCNVs(chldrnIds,minGn):
    PSD = {pId:0 for pId in chldrnIds}
    for e in CNVs:
        if len(e.gns)<minGn: 
            continue
        for pid in e.pids:
            if pid in PSD:
                PSD[pid] += 1
    return {pid:(c,1) for pid,c in PSD.items()}

clf()
geneNumbers = [0,1,2,3,4,5]
xlbls = []
for gi,gn in enumerate(geneNumbers):
    prbC = countCNVs(prbs,gn)
    sibC = countCNVs(sibs,gn)
    sReal, sBtstrp = compareN(prbC,sibC,bootstrapI=1000)
    bcIR = empStats(sReal, sBtstrp, 'IR')
    plot(gn,sReal.IR,'*r')
    plot([gn,gn],[bcIR.left95,bcIR.right95],'r')
    xlbls.append("%d\n(%d, %d)" % (gn, sReal.Na,sReal.Nu))
xticks(geneNumbers,xlbls)
ylabel('percent contributory')
xlabel('min CNV gene number\n(number of CNVs in affected and unaffected children)')
gcf().set_size_inches(6,4)
tight_layout()
if outFigName:
    gcf().savefig(outFigName)
else:
    show()
