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


IFN = "intronic_result_table.txt"
if len(sys.argv) > 2:
    IFN = sys.argv[2] 

RS = genfromtxt(IFN, delimiter='\t',dtype=None,names=True, case_sensitive=True)

RSD = {}
for R in RS:
    k = (R['set'],R['eventType'],R['intronType'])
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
        affN = unwR['sRealNa']
        affX = unwR['sRealXa']
        unaN = unwR['sRealNu']
        unaX = unwR['sRealXu']

        dlt = unwR['sRealdelta']
        expAffN = unwR['sRealENa']

        AD = unwR['sRealAD']
        ADl = unwR['bcADleft95']
        ADr = unwR['bcADright95']

        IR = unwR['sRealIR']
        IRl = unwR['bcIRleft95']
        IRr = unwR['bcIRright95']

        nGenesS = "{:,}".format(int(unwR['setGeneNumber']))
        cs = [stT,nGenesS,etpOut] + \
            ["{:,}".format(int(v)) for v in unaX,unaN,affX,affN] + \
            ["{:,.1f}".format(expAffN), 
             "%.1f" % dlt, \
             pV2Str(float(unwR['esADpvOne'])), \
             "%.2f%% (%.2f-%.2f)" % (AD,ADl,ADr), \
             "%.1f%% (%.1f-%.1f)" % (IR,IRl,IRr) \
            ]

        print "\t".join(map(str,cs))
