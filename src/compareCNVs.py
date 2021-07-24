#!/usr/bin/env python

from diData import loadEVS, resDir
from collections import Counter

CNVs = loadEVS(['De novo CNV in SSC and AGRE'])

A = {(e.pid,e.variant,e.location) for e in CNVs }

F = open(resDir + "/Supplementary Data 6. CNV details.txt")
hl = F.readline()
print hl,
hcs = hl.strip("\n\r").split("\t")
inSCI = hcs.index('in_summary')

B = set([])
for l in F:
    cs = l.strip("\n\r").split("\t")
    assert len(cs) == len(hcs)
    cnv = dict(zip(hcs,cs))
    pid = cnv['sampleID']

    if cnv['if_merged'] == "FALSE" or \
       cnv['parents_dels'] == 'NA' or cnv['parents_dups'] == 'NA' or int(cnv['parents_dels']) > 5 or int(cnv['parents_dups']) > 5 or \
       cnv['ploidy_measure'] == 'NA' or float(cnv['ploidy_measure']) < 0.85 or float(cnv['ploidy_measure']) > 1.15 or \
       pid in outlierPersons:
        cs[inSCI] = 'FALSE' 
    else:
        cs[inSCI] = 'TRUE' 
        B.add((pid,cnv['polarity'],"chr%s:%s-%s" % (cnv['chr'],cnv['start'],cnv['end'])))
    print "\t".join(cs)
F.close()
