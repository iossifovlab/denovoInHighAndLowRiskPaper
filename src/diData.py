#!/usr/bin/env python

import sys,glob,os
from pylab import *
from collections import defaultdict

resDir = "./input"
if 'DI_DATA_RES_DIR' in os.environ:
    resDir = os.environ['DI_DATA_RES_DIR']

eventFiles = glob.glob('%s/*_denovo_table.txt' % (resDir))
if 'DI_DATA_EVENT_FILES' in os.environ:
    eventFiles = os.environ['DI_DATA_EVENT_FILES'].split(":")

print "resDir:", resDir
print "eventFiles:", eventFiles

try:
    DT
except NameError:
    DT = genfromtxt(resDir + "/gene_table.txt", delimiter='\t',dtype=None,names=True, case_sensitive=True,encoding=None)

class Child:
    pass
class Family:
    pass

try:
    CHILDREN
except NameError:
    CHILDREN = genfromtxt(resDir + "/children_table.txt", delimiter='\t',dtype=None,names=True, case_sensitive=True,encoding=None)

persons = defaultdict(int)
outlierPersons = defaultdict(int)
families = defaultdict(int)
afs2Role = {"affected":"prb", "unaffected":"sib"}
for R in CHILDREN:
    fId = str(R['familyId'])
    p = Child()
    p.fId = fId
    p.pId,p.affectedStatus,p.gndr,p.coll,p.outlier,p.dnaSource = [R[a] for a in "personId,affected_status,sex,collection,cell_line_drift_outlier,dna_source".split(",")]
    p.role = afs2Role[p.affectedStatus]
    p.atts = R

    assert p.pId not in persons
    assert p.pId not in outlierPersons 
    if p.outlier == 1:
        outlierPersons[p.pId] = p
    else:
        persons[p.pId] = p

    if fId not in families:
        families[fId] = Family()
        families[fId].fId = fId
        families[fId].chldrn = []
        families[fId].outlierChldrn = []
    if p.outlier == 1:
        families[fId].outlierChldrn.append(p)
    else:
        families[fId].chldrn.append(p)

for fd in families.values():
    fd.coll ,= set([p.coll for p in fd.chldrn + fd.outlierChldrn])

class Denovo:
    pass

try:
    EVS
except NameError:
    personSet = set(persons)
    EVS = [] 
    for fffnnn in eventFiles:
        print >>sys.stderr, "loading", fffnnn, "..."
        try:
            EVD = genfromtxt(fffnnn, delimiter='\t',dtype=None,names=True, case_sensitive=True,encoding=None)    
        except Exception:
            print >>sys.stderr, "FAILED TO LOAD", fffnnn
            continue
        for r in EVD:
            dn = Denovo()
            dn.fmid =  str(r['familyId'])

            dn.pids = r['personIds'].split(",")
            assert set(dn.pids).issubset(personSet)
            if len(dn.pids)==1: dn.pid ,= dn.pids

            dn.coll ,= set([persons[pid].coll for pid in dn.pids])

            dn.gns = set([])
            if r['genes']:
                dn.gns = set(r['genes'].split(","))

            if 'variant_type' in r.dtype.fields:
                dn.vtype = r['variant_type']
            else:
                dn.vtype = 'CNV'

            dn.variant = r['variant']
            dn.location = r['location']
           
            if 'effect_type' in r.dtype.fields:
                dn.eff = r['effect_type']
            else:
                dn.eff = dn.variant

            dn.genomicRegion = r['genomic_region']

            dn.atts = r
            
            EVS.append(dn)

if __name__ == "__main__":
    print 'hi'

