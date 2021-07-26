#!/usr/bin/env python

from pylab import *
from pV2Str import pV2Str

IFN = "resTab-smallScale.txt"
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
    k = (R['group'],R['backgroundGroup'],R['effect'])
    RSD[k] = R

hcs1 = ['','SSC unaffected','','affected']
hcs2 = ['group','synonymous number', 'LGD number', 'synonymous number', 'LGD number', 'expected LGD number', 'delta', 'pval', 'AD', 'PC']
print "\t".join(hcs1)
print "\t".join(hcs2)
for grp in ['SSC affected', 'AGRE affected']:
    R = RSD[grp,'SSC unaffected','LGD']
    cs = [grp] + \
         ['{:,}'.format(int(R[a])) for a in ['Nu', 'Su', 'Na', 'Sa']] + \
         ["%.1f" % (float(R[a])) for a in ['ESa','delta']] + \
         [pV2Str(float(R['AD.pvOneAn']))] + \
         ["%.2f%% (%.2f-%.2f)" % tuple([float(R[a]) for a in ['AD','AD.left95','AD.right95']])] + \
         ["%.1f%% (%.1f-%.1f)" % tuple([float(R[a]) for a in ['PC','PC.left95','PC.right95']])]
    print "\t".join(cs)


