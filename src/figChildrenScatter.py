#!/usr/bin/env python

from pylab import *
from twodnorm import *
from params import *
import matplotlib.gridspec as gridspec
from matplotlib.transforms import Bbox
from scipy import stats

CH = genfromtxt("children_table.txt", delimiter='\t',dtype=None,names=True, case_sensitive=True,encoding=None)

xs = CH['nDNSnvs']/CH['power']
ys = CH['meanDNSnvsAltAlleleRatio']

xswb = xs[CH['dnaSource']=='WB']
yswb = ys[CH['dnaSource']=='WB']

XWB = vstack([xswb,yswb]).T
mns = XWB.mean(axis=0)
cvM = cov(XWB.T)

outlierConf = float(params['children.binormal.conf.cutoff'])

XALL = vstack([xs,ys]).T
# cfs = assignConf(mns,cvM,XALL)
cfs = CH['biVariantNormConf']


def xTrns(x):
    x = array(x)
    # return x
    C = 150 
    x[x>C] = C + (x[x>C]-C)/10.
    x[x>C] = C + C*log(x[x>C]-C)/log(1900)
    return x

clf()
GSS = 5 
gs = gridspec.GridSpec(GSS, GSS)
ax_main = subplot(gs[1:GSS, :(GSS-1)])

cs = []
for R in CH:
    c = 'g'
    if R['collection'] == 'AGRE':
        if R['dnaSource'] == 'WB':
            c = 'r'
        else:
            c = 'b'
    cs.append(c)
cs = array(cs)
        

# gca().add_patch(confEllipse(mns,cvM,0.05,ec='m',lw=1))
# gca().add_patch(confEllipse(mns,cvM,0.01,ec='g',lw=1))
gca().add_patch(confEllipse(mns,cvM,outlierConf,ec='k',lw=1))

# ss = 20*((cfs>outlierConf).astype(int)) + 1
ss = 5*ones(len(cfs))

xxx = XALL[:,0]
scatter(xTrns(xxx),XALL[:,1],c=cs,s=ss)
xtcks = [0,50,100,150,200,500,2000]
xticks(xTrns(xtcks),map(str,xtcks))

SWI = CH['collection']=='SSC'
AWI = logical_and(CH['collection']=='AGRE',CH['dnaSource']=='WB')
ALI = logical_and(CH['collection']=='AGRE',CH['dnaSource']=='LCL')

OKI = CH['biVariantNormConf'] > outlierConf

def pltD(a,c,f,t,onY=False,):
    density = stats.kde.gaussian_kde(a)
    x = linspace(f,t,1000)
    if onY:
        plot(density(x), x,c)
    else:
        plot(x, density(x),c)
plot([150,150],ylim(),'k:')
xlabel('power adjusted number of de novo substitutions')
ylabel('mean alternative allele ratio')

# Hide the right and top spines
gca().spines['right'].set_visible(False)
gca().spines['top'].set_visible(False)


ax_xDist = subplot(gs[0, :(GSS-1)],sharex=ax_main)
xxxx = xTrns(xxx)
M = xTrns(4000)
pltD(xxxx[SWI],'g',0,M)
pltD(xxxx[AWI],'r',0,M)
pltD(xxxx[ALI],'b',0,M)
axis('off')

legend([  'SSC WB %d/%d'  % (logical_and(SWI,OKI).sum(),  SWI.sum()),
         'AGRE WB %d/%d'  % (logical_and(AWI,OKI).sum(),  AWI.sum()), 
         'AGRE LCL %d/%d' % (logical_and(ALI,OKI).sum(),  ALI.sum())])

ax_yDist = subplot(gs[1:GSS, GSS-1],sharey=ax_main)
yyy = XALL[:,1]
pltD(yyy[SWI],'g',0.27,0.55,True)
pltD(yyy[AWI],'r',0.27,0.55,True)
pltD(yyy[ALI],'b',0.27,0.55,True)
axis('off')

gcf().set_size_inches(8,8)
tight_layout()

bp = ax_yDist.get_position().get_points()
print "y", bp
bp[0,0] = 0.75
ax_yDist.set_position(Bbox(bp))

bp = ax_xDist.get_position().get_points()
print "x", bp
bp[0,1] = 0.75
ax_xDist.set_position(Bbox(bp))

gcf().savefig("children_scatter.png")
