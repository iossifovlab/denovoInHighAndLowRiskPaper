#!/usr/bin/env python

from numpy import sqrt, arctan, pi, argsort
from scipy.stats import multivariate_normal as MVN,chi2
import matplotlib.patches as mpatches
from numpy.linalg import eig

def sortedEig(cv):
    l,uv = eig(cv)
    # sort by eigenvalue 
    srti = argsort(l)[::-1]
    l = l[srti]
    uv = uv[:,srti]
    return l,uv

def confEllipse(mn,cv,conf,fill=False,ec='k',fc='b',lw=None):
    s = chi2.isf(conf,2)

    l,uv = sortedEig(cv)

    # elipse 
    W,H = 2*sqrt(s*l)
    a = 180 * arctan (uv[1,0]/uv[0,0]) / pi
    return mpatches.Ellipse(mn,W,H,a,fill=False,ec=ec,fc=fc,lw=lw)

def assignConf(mn,cv,R):
    l,uv = sortedEig(cv)

    RT = (R-mn).dot(uv)
    v = ((RT**2 / l)).sum(axis=1)
    return chi2.sf(v,2)

# return (firstXs,firstYs),(secondXs,secondYs)
def getAxesVectors(mn,cv):
    l,uv = sortedEig(cv)
    axp = uv * l + mn[np.newaxis].T
    return ( [mn[0], axp[0,0]], [mn[1],axp[1,0]]), \
           ( [mn[0], axp[0,1]], [mn[1],axp[1,1]]) 

if __name__ == "__main__":
    from pylab import *
    mns = array([15.,30.])
    cv = 4.0
    cvM = array([[4.,cv],[cv,9.]])

    R = MVN.rvs(mns,cvM,size=10000)

    clf()
    plot(R[:,0],R[:,1],'.b')
    gca().add_patch(confEllipse(mns,cvM,0.05,ec='k',lw=3))
    gca().add_patch(confEllipse(mns,cvM,0.01,ec='g',lw=3))
    gca().add_patch(confEllipse(mns,cvM,0.001,ec='m',lw=3))
    cfs = assignConf(mns,cvM,R)
    pltConf = 0.05
    plot(R[cfs<pltConf,0],R[cfs<pltConf,1],'.r')
    
    (fXs,fYs),(sXs,sYs) = getAxesVectors(mns,cvM)
    plot(fXs,fYs,'r',lw=5)
    plot(sXs,sYs,'g',lw=5)

    show()

