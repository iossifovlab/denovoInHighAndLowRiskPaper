#!/usr/bin/env python

from pylab import *
from pV2Str import pV2Str

mtSets = [
    ('Y',"all genes"),
    ('Y',"autism LGD"),
    ('N',"all NDD LGD"),
    ('N',"autism missense"),
    ('N',"autism synonymous"),
]

itp = 'inter-coding_intronic'

if len(sys.argv) > 1:
    itp = sys.argv[1] 


IFN = "resTab-intronicPeripheral.txt"
if len(sys.argv) > 2:
    IFN = sys.argv[2] 

RS = genfromtxt(IFN, delimiter='\t',dtype=None,names=True, case_sensitive=True,encoding=None)

RSD = {}
for R in RS:
    k = (R['geneSet'],R['eventType'],R['regionType'])
    RSD[k] = R

hcs1 = ['','','','SSC unaffected', '', 'SSC affected']
hcs2 = ",gene,event,normalization,functional,normalization,functional,expected functional,,,,".split(",")
hcs3 = "set,number,type,number,number,number,number,number,delta,pvalue,AD,PC".split(",")

print "\t".join(hcs1)
print "\t".join(hcs2)
print "\t".join(hcs3)

for etp,etpOut in zip(['indel', 'sub'],['IID','ISB']):
    for stSub,stK in mtSets:
        if etp == 'sub' and stSub != 'Y':
            continue

        if stK == "all genes":
            stT = stK
        else:
            stT = stK + " targets"
        unwR = RSD[(stK,etp,itp)]
        affN = unwR['Sa']
        affX = unwR['Na']
        unaN = unwR['Su']
        unaX = unwR['Nu']

        dlt = unwR['delta']
        expAffN = unwR['ESa']

        AD = unwR['AD']
        ADl = unwR['ADleft95']
        ADr = unwR['ADright95']

        PC = unwR['PC']
        PCl = unwR['PCleft95']
        PCr = unwR['PCright95']

        nGenesS = "{:,}".format(int(unwR['geneSetGeneNumber']))
        cs = [stT,nGenesS,etpOut] + \
            ["{:,}".format(int(v)) for v in unaX,unaN,affX,affN] + \
            ["{:,.1f}".format(expAffN), 
             "%.1f" % dlt, \
             pV2Str(float(unwR['ADpvOne'])), \
             "%.2f%% (%.2f-%.2f)" % (AD,ADl,ADr), \
             "%.1f%% (%.1f-%.1f)" % (PC,PCl,PCr) \
            ]

        print "\t".join(map(str,cs))
