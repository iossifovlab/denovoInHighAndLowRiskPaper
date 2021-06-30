#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')

from pylab import *
from scipy.stats import beta


CDS = '''size	ssc_aut	ssc_sib	agre_aut
1 kbp	8738	8689	8079
2 kbp	9611	9617	9493
3 kbp	9721	9703	9642
4 kbp	9752	9742	9725
5 kbp	9750	9731	9736
6 kbp	9716	9723	9752
7 kbp	9768	9756	9775
8 kbp	9763	9747	9738
9 kbp	9766	9737	9780
10 kbp	9801	9782	9754'''

'''
20 kbp	9776	9744	9790
30 kbp	9797	9792	9761
40 kbp	9803	9807	9803
50 kbp	9797	9765	9821
60 kbp	9802	9802	9801
70 kbp	9810	9795	9814
80 kbp	9831	9805	9789
90 kbp	9817	9806	9799
100 kbp	9825	9811	9796'''



CDS_old = '''size	agre	ssc
2kbp	8806	9117
3kbp	9354	9506
4kbp	9641	9636
5kbp	9726	9719
6kbp	9742	9738
7kbp	9760	9752
8kbp	9778	9781
9kbp	9791	9770
10kbp	9799	9768
20kbp	9817	9798
30kbp	9820	9802
40kbp	9844	9805
50kbp	9845	9840
60kbp	9839	9825
70kbp	9854	9846
80kbp	9864	9831
90kbp	9853	9841
100kbp	9874	9842'''

ls = CDS.split("\n");
hcs = ls[0].split("\t")
    
def plotFCI(n,x,c):
    bd = beta(n,10000-n)
    lid = plot([x,x],[bd.ppf(0.05),bd.isf(0.05)],c,lw=2,alpha=0.5)
    plot([x],n/10000.,c + '*')
    return lid[0]

clf()
xls = []
for x,l in enumerate(ls[1:]):
    cs = l.split('\t')
    assert len(cs) == len(hcs)

# CDS = '''size	ssc_aut	ssc_sib	agre_aut

    cnvSize,sscAutOK,sscSibOK,agreAutOK = cs
    agreAutOK = int(agreAutOK)
    sscAutOK = int(sscAutOK)
    sscSibOK = int(sscSibOK)
    aal = plotFCI(agreAutOK,x,'k')
    sal = plotFCI(sscAutOK,x,'r')
    ssl = plotFCI(sscSibOK,x,'g')
    xls.append(cnvSize)
# xticks(arange(len(xls)),xls,fontsize=8,rotation=-90)
xticks(arange(len(xls)),xls)
yl = ylim()
ylim([yl[0],1])
grid(b=True,axis='x')
xlabel('CNV size')
ylabel('power')
legend([sal,ssl,aal],['SSC affected','SSC unaffected','AGRE affected'])
gcf().set_size_inches(7,3)
tight_layout()
if len(sys.argv) > 1:
    gcf().savefig(sys.argv[1])
else:
    show()
