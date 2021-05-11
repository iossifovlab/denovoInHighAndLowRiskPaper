#!/usr/bin/env python

from pylab import *
import matplotlib.gridspec as gridspec
from matplotlib.transforms import Bbox
from scipy import stats
import twodnorm
from diData import *

CH = persons.values() + outlierPersons.values()
NCH = len(CH)

xs = array([pd.atts[CN('number of dn snvs')]/pd.atts[CN('power to detect dn snvs')] for pd in CH])
ys = array([pd.atts[CN('mean dn snvs alt allele ratio')] for pd in CH])

SWI = array([pd.coll == 'SSC' for pd in CH])
AWI = array([pd.coll == 'AGRE' and pd.dnaSource == 'WB' for pd in CH])
ALI = array([pd.coll == 'AGRE' and pd.dnaSource == 'LCL' for pd in CH])

WBI = array([pd.dnaSource == 'WB' for pd in CH])

OKI = array([pd.outlier != 1 for pd in CH])

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


##        
## draw the 2-d confidence 0.001 ellipse
##
outlierConf = 0.001

# step 1. fit the 2d normal distribution to the WB samples
xswb = xs[WBI]
yswb = ys[WBI]
XWB = vstack([xswb,yswb]).T
mns = XWB.mean(axis=0)
cvM = cov(XWB.T)

# step 2. plot the ellipse
gca().add_patch(twodnorm.confEllipse(mns,cvM,outlierConf,ec='k',lw=1))

# verity that the outlier samples are the ones with confidence level less than 0.001
XALL = vstack([xs,ys]).T
cfs = twodnorm.assignConf(mns,cvM,XALL)
assert ( OKI == (twodnorm.assignConf(mns,cvM,vstack([xs,ys]).T) > outlierConf) ).all()

##########################
##########################
##########################

ss = 5*ones(NCH)
# ss[OKI] = 20

cs = array(['g']*NCH)
cs[AWI] = 'r'
cs[ALI] = 'b'

scatter(xTrns(xs),ys,c=cs,s=ss)
xtcks = [0,50,100,150,200,500,2000]
xticks(xTrns(xtcks),map(str,xtcks))

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
xxxx = xTrns(xs)
M = xTrns(4000)
pltD(xxxx[SWI],'g',0,M)
pltD(xxxx[AWI],'r',0,M)
pltD(xxxx[ALI],'b',0,M)
axis('off')

legend([  'SSC WB %d/%d'  % (logical_and(SWI,OKI).sum(),  SWI.sum()),
         'AGRE WB %d/%d'  % (logical_and(AWI,OKI).sum(),  AWI.sum()), 
         'AGRE LCL %d/%d' % (logical_and(ALI,OKI).sum(),  ALI.sum())])

ax_yDist = subplot(gs[1:GSS, GSS-1],sharey=ax_main)
pltD(ys[SWI],'g',0.27,0.55,True)
pltD(ys[AWI],'r',0.27,0.55,True)
pltD(ys[ALI],'b',0.27,0.55,True)
axis('off')

gcf().set_size_inches(8,8)
tight_layout()

bp = ax_yDist.get_position().get_points()
bp[0,0] = 0.75
ax_yDist.set_position(Bbox(bp))

bp = ax_xDist.get_position().get_points()
bp[0,1] = 0.75
ax_xDist.set_position(Bbox(bp))

if len(sys.argv) > 1:
    gcf().savefig(sys.argv[1])
else:
    show()
