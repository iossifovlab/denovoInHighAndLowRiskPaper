#!/usr/bin/env python

from diData import resDir, persons, outlierPersons
from collections import Counter

F = open(resDir + "/Supplementary Data 6. CNV details.txt")
hl = F.readline()
print hl,
hcs = hl.strip("\n\r").split("\t")
assert not set(persons) & set(outlierPersons)
knownPersons = set(persons) | set(outlierPersons)

for l in F:
    cs = l.strip("\n\r").split("\t")
    assert len(cs) == len(hcs)
    R = dict(zip(hcs,cs))
    pid = R['sampleID']

    if pid not in knownPersons:
        continue
    print l, 
F.close()
