#!/usr/bin/env python

from pylab import *
from diData import *
from scipy import stats
import sys

def pltD(aAll,c,f,t,lw=1,onY=False,):
    a = array(aAll)[logical_not(isnan(aAll))] 
    density = stats.kde.gaussian_kde(a)
    x = linspace(f,t,1000)
    if onY:
        plot(density(x), x,c,lw=lw)
    else:
        plot(x, density(x),c,lw=lw)
    return len(a)

chldSetDef = [['SSC','affected','r'],['SSC','unaffected','g'], ['AGRE','affected', 'k']]

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
