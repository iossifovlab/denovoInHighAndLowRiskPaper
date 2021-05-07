#!/usr/bin/env python

from pylab import *

   

IFN = "merged_simplex_multiplex_results.txt"
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
    k = (R['foreground group'],R['background group'])
    RSD[k] = R

def pV2Str(pv):
    if pv < 0.0001:
        return "%.0e" % (pv)
    elif pv < 0.001:
        return "%.5f" % (pv)
    elif pv < 0.01:
        return "%.4f" % (pv)
    elif pv < 0.1:
        return "%.3f" % (pv)
    else:
        return "%.2f" % (pv)

hcs = ['group','events number', 'expected events number', 'delta', 'pval', 'AD']
print "\t".join(hcs)
for grp in ['SSC affected', 'AGRE affected']:
    R = RSD[grp,'SSC unaffected']
    cs = [grp, R['B']]

    cs += ["%.1f" % (float(R[a])) for a in ['EB','delta']] + \
          [pV2Str(float(R['AD.pvOneAn']))] + \
          ["%.2f%% (%.2f-%.2f)" % tuple([float(R[a]) for a in ['AD','AD.left','AD.right']])] 
    print "\t".join(cs)


