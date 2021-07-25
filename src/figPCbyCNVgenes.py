#!/usr/bin/env python

import sys
import matplotlib as mpl
if len(sys.argv) > 1:
    mpl.use('Agg')

from pylab import *
from diData import persons, loadEVS
from methods import compare_subject_variant_class_in_two_groups_using_normalization_variant_class, empStats

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
    sReal, sBtstrp = compare_subject_variant_class_in_two_groups_using_normalization_variant_class(prbC,sibC,bootstrapI=1000)
    bcPC = empStats(sReal, sBtstrp, 'PC')
    plot(gn,sReal.PC,'*r')
    plot([gn,gn],[bcPC.left95,bcPC.right95],'r')
    xlbls.append("%d\n(%d, %d)" % (gn, sReal.Sa,sReal.Su))
xticks(geneNumbers,xlbls)
ylabel('percent contributory')
xlabel('min CNV gene number\n(number of CNVs in affected and unaffected children)')
gcf().set_size_inches(6,4)
tight_layout()
if outFigName:
    gcf().savefig(outFigName)
else:
    show()
