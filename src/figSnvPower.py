#!/usr/bin/env python

import sys
import matplotlib as mpl
if len(sys.argv) > 1:
    mpl.use('Agg')

from pylab import *
from diData import persons, CN
from scipy import stats
from collections import defaultdict

def pltD(aAll,c,f,t,lw=1,onY=False,):
    a = array(aAll)[logical_not(isnan(aAll))] 
    density = stats.kde.gaussian_kde(a)
    x = linspace(f,t,1000)
    if onY:
        plot(density(x), x,c=c,lw=lw)
    else:
        plot(x, density(x),c=c,lw=lw)
    return len(a)

chldSetDef = [
    ['SSC' , 'affected'  , [213/255., 94 /255., 0  /255.]],
    ['SSC' , 'unaffected', [0  /255., 158/255., 115/255.]], 
    ['AGRE', 'affected'  , [0  /255., 0  /255., 0  /255.]]
]

figure()
GRPP = defaultdict(list)
for pd in persons.values():
    GRPP[pd.coll,pd.affectedStatus].append(pd.atts[CN('power to detect dn snvs')])

lgnds = []
for coll,affSt,c in chldSetDef:
    pwrs = GRPP[coll,affSt]
    pltD(pwrs,c,0.4,1.0)
    lgnds.append("%s %s %d" % (coll,affSt,len(pwrs)))
legend(lgnds,loc='upper left')
xlabel('power')
ylabel('pdf')
yticks([])

gcf().set_size_inches(5,5)
tight_layout()
if len(sys.argv) > 1:
    gcf().savefig(sys.argv[1])
else:
    show()
