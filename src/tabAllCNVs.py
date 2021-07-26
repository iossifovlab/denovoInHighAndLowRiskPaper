#!/usr/bin/env python

from pylab import *
from pV2Str import pV2Str

IFN = "resTab-CNVs.txt"
if len(sys.argv) > 1:
    IFN = sys.argv[1] 

RS = []
IF = open(IFN)
hcs = IF.readline().strip("\n\r").split("\t")
for l in IF:
    if l.startswith("#"): continue
    l = l.strip("\n\r")
    if not l: continue
    cs = l.split("\t")
    assert len(cs) == len(hcs)
    R = dict(zip(hcs,cs))
    RS.append(R)
IF.close()

RSD = {}
for R in RS:
    k = (R['section'],R['effect'],R['variant'],R['group'],R['backgroundGroup'])
    RSD[k] = R

hcs1 = ['','','SSC unaffected','','affected']
hcs2 = ['effect', 'group','CNV number', 'CNV rate', 'CNV number', 'CNV rate', 'expected CNVs number', 'delta', 'pval', 'AD', 'PC']
print "\t".join(hcs1)
print "\t".join(hcs2)
# for vTp in ['cnvs','deletion','duplication']:
for effT in ['all','coding','intergenic','genic noncoding']:
    for grp in ['SSC affected', 'AGRE affected']:
        R = RSD['ALL EVENTS',effT,'cnvs',grp,'SSC unaffected']
        cs = [effT,grp]
        for s in 'ua':
            N = int(R['S' + s])
            C = int(R['N' + s])
            cs += [str(N), '%.3f' % (float(N)/C)]

        cs += ["%.1f" % (float(R[a])) for a in ['ESa','delta']] + \
              [pV2Str(float(R['AD.pvOneAn']))] + \
              ["%.2f%% (%.2f-%.2f)" % tuple([float(R[a]) for a in ['AD','AD.left95','AD.right95']])] + \
              ["%.1f%% (%.1f-%.1f)" % tuple([float(R[a]) for a in ['PC','PC.left95','PC.right95']])]
        print "\t".join(cs)


