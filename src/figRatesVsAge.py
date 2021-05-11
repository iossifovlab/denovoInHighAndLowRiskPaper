#!/usr/bin/env python

from pylab import *
from diData import *
from scipy import stats
from collections import defaultdict

EVS = loadEVS(['SSC_small_denovo'])

CNT = defaultdict(lambda : {pid:0 for pid,pd in persons.items() if pd.coll == 'SSC'})
 
for e in EVS:
    parent = e.atts[CN('from parent')]
    if not parent: continue
    if e.location.startswith('chrX'): continue

    tp = "sub" if e.vtype == "sub" else "indel"

    for pid in e.pids:
        CNT[tp,parent][pid] += 1


def p(tp,parent):
    nns = []
    ags = []

    for pid,n in CNT[tp,parent].items():
        ag = persons[pid].atts[CN('father age at birth' if parent == 'dad' else 'mother age at birth')]
        if isnan(ag): continue
        if ag < 15. or ag > 60: continue
        ags.append(ag)
        nns.append(n/persons[pid].atts[CN('power to detect dn snvs')])

    ags = array(ags)
    nns = array(nns)

    slope, intercept, rvalue, pvalue, stderr = stats.linregress(ags,nns)
    plot(ags,nns,'.')
    mx = array([15.0,60.0])
    my = slope * mx + intercept
    plot(mx,my)
    pvt = "%.0e" % (pvalue)
    g,h = pvt.split("e")
    hi = int(h)
    if hi > -6.:
        tpvt = "%s*$10^{-%s}$" % (g,h.lstrip('-0'))
    else:
        tpvt = "$10^{-%s}$" % (h.lstrip('-0'))
    print pvt

    text(15.,0.85*ylim()[1],'slope: %.1f / 10 years\np-value: %s' % (slope * 10.,tpvt))

clf()
subplot(2,2,1)
ylim([0,50])
p('sub','dad')
ylabel('power adjusted number of\nphased dn substitutions')
title('Father')

subplot(2,2,2)
ylim([0,50])
p('sub','mom')
title('Mother')

subplot(2,2,3)
ylim([0,12])
p('indel','dad')
ylabel('power adjusted number of\nphased dn indels')
xlabel("father's age at birth")

subplot(2,2,4)
ylim([0,10])
p('indel','mom')
xlabel("mother's age at birth")

gcf().set_size_inches(6.5,6.5)
tight_layout()
gcf().savefig("figRatesVsAge.png")
# show()
