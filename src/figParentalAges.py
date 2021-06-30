#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')

from diData import *
from scipy import stats

GRPSN = defaultdict(list)
GRPMA = defaultdict(list)
GRPDA = defaultdict(list)

for pd in persons.values():
    k = (pd.coll,pd.affectedStatus)

    GRPSN[k].append(pd.atts[CN('number of dn snvs')]/pd.atts[CN('power to detect dn snvs')])
    GRPMA[k].append(pd.atts[CN('mother age at birth')])
    GRPDA[k].append(pd.atts[CN('father age at birth')])

def pltD(aAll,c,f,t,lw=1,onY=False,):
    a = array(aAll)[logical_not(isnan(aAll))] 
    density = stats.kde.gaussian_kde(a)
    x = linspace(f,t,1000)
    if onY:
        plot(density(x), x,c,lw=lw)
    else:
        plot(x, density(x),c,lw=lw)
    return len(a)


clf()
spi = 0
chldSetDef = [['SSC','affected','r'],['SSC','unaffected','g'], ['AGRE','affected', 'k']]
for par,GS in zip(['father','mother'],[GRPDA,GRPMA]):
    spi += 1
    subplot(2,2,spi)
    lgnds = []
    for coll,affSt,c in chldSetDef:
        agesAll = array(GS[coll,affSt])
        ages = agesAll[logical_not(isnan(agesAll))] 
        pltD(ages,c,10,100)
        lgnds.append('%s %s %d/%d' % (coll, affSt, len(ages),len(agesAll)))
    title(par + "'s age at birth")
    legend(lgnds,fontsize=8,loc='upper right')
    # ylim([-0.1,0.2])
    yticks([])
    xlim([10,75])

    spi += 1
    subplot(2,2,spi)
    lgnds = []
    for ci,(coll,affSt,c) in enumerate(chldSetDef):
        agesAll = array(GS[coll,affSt])
        ages = agesAll[logical_not(isnan(agesAll))] 
        mn = mean(ages)
        se = 2*stats.sem(ages)
        plot([ci],[mn],c+'o')
        plot([ci,ci],[mn-se,mn+se],c+'-')
        plot([0,2],[mn,mn],c+':')
        lgnds.append('%s %s' % (coll, affSt))
    xticks([0,1,2],lgnds)
    ylabel('age')
    title(par + "'s mean age at birth")
gcf().set_size_inches(8,6)
tight_layout()
if len(sys.argv) > 1:
    gcf().savefig(sys.argv[1])
else:
    show()
