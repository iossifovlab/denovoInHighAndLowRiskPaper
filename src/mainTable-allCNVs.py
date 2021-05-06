#!/usr/bin/env python

from pylab import *

   

IFN = "CNV_result_table.txt"
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

hcs1 = ['','','SSC unaffected','','affected']
hcs2 = ['effect', 'group','CNV number', 'CNV rate', 'CNV number', 'CNV rate', 'expected CNVs number', 'delta', 'pval', 'AD', 'PC']
print "\t".join(hcs1)
print "\t".join(hcs2)
# for vTp in ['cnvs','deletion','duplication']:
for effT in ['all','coding','intergenic','genic other']:
    for grp in ['SSC affected', 'AGRE affected']:
        R = RSD['ALL EVENTS',effT,'cnvs',grp,'SSC unaffected']
        cs = [effT,grp]
        for s in 'ua':
            N = int(R['sReal.N' + s])
            C = int(R['sReal.X' + s])
            cs += [str(N), '%.3f' % (float(N)/C)]

        cs += ["%.1f" % (float(R[a])) for a in ['sReal.ENa','sReal.delta']] + \
              [pV2Str(float(R['esAD.pvOneAn']))] + \
              ["%.2f%% (%.2f-%.2f)" % tuple([float(R[a]) for a in ['sReal.AD','bcAD.left95','bcAD.right95']])] + \
              ["%.1f%% (%.1f-%.1f)" % tuple([float(R[a]) for a in ['sReal.IR','bcIR.left95','bcIR.right95']])]
        print "\t".join(cs)


