from numpy import array, sort, nan
from numpy.random import permutation, randint
import scipy.stats as st

def compareOneSubjectClass(countsInAffected,countsInUnaffected,nullI=0,bootstrapI=0):
    '''Compares a subject variant class in two groups of children 
       using a normalization variant class. The two groups of children
       are referred to as affected and unaffected. 
       The countsInAffected and the countsInUnaffected parameters have the form:
        {
            "personId1":["number of subject class variants", "number of normalization class variants],
            "personId2":["number of subject class variants", "number of normalization class variants],
            ....
        }

        Returns an instance of the VStats class that contains statistics of comparison of the two groups.
        See below fo details and the manuscript to details of the fields. If the nullI > 0, the function will
        also return an array of VStats for nullI random permutation of the affected statuses of the children.
        Similarly, if the bootstrapI > 0, the function will also return an array of VStats for 
        bootstrapI permutation iteration.
    '''
    SV = []
    NV = []
    CV = []
    for ci,c in enumerate([countsInUnaffected,countsInAffected]):
        for pid,(n,x) in sorted(c.items()):
            SV.append(n) 
            NV.append(x) 
            CV.append(ci) 
    SV = array(SV)
    NV = array(NV)
    CV = array(CV)

    def vStats(SV,NV,CV):
        class VStats:
            pass
        s = VStats() 

        s.Ca = (CV==1).sum()
        s.Cu = (CV==0).sum()
        s.Sa = SV[CV==1].sum()
        s.Su = SV[CV==0].sum() 
        s.Na = NV[CV==1].sum()
        s.Nu = NV[CV==0].sum()

        s.SNa = float(s.Sa) / s.Na
        s.SNu = float(s.Su) / s.Nu

        s.ESa = s.SNu * s.Na
        s.delta = s.Sa - s.ESa

        s.AD = 100. * s.delta / s.Ca 

        if s.Sa:
            s.PC = 100. * s.delta / s.Sa 
        else:
            s.PC = 100. 

        s.RSa = float(s.Sa)/s.Ca
        s.RSu = float(s.Su)/s.Cu
        s.RNa = float(s.Na)/s.Ca
        s.RNu = float(s.Nu)/s.Cu

        return s

    sReal = vStats(SV,NV,CV)
    r = [sReal]
    if nullI:
        r.append([vStats(SV,NV,permutation(CV)) for e in xrange(nullI)])
    if bootstrapI:
        bts = []
        for e in xrange(bootstrapI):
            ris = randint(len(SV),size=len(SV))
            bts.append(vStats(SV[ris],NV[ris],CV[ris]))
        r.append(bts)
    return r 


def compareTwoSubjectClasses(affectedCounts1,unaffectedCounts1,affectedCounts2,unaffectedCounts2,nullI=0,bootstrapI=0):
    '''Compares joinly two subject variant classes (1 and 2) in two groups of 
       children (affected and unaffected) using 
       separate normalization variant classes.

       The arguments and the return values are analogous to the function compareOneSubjectClass above.
    '''
    assert affectedCounts1.keys() == affectedCounts2.keys()
    assert unaffectedCounts1.keys() == unaffectedCounts2.keys()

    S1V = []
    N1V = []
    S2V = []
    N2V = []
    CV = []
    for ci,(c1,c2) in enumerate([(unaffectedCounts1,unaffectedCounts2),(affectedCounts1,affectedCounts2)]):
        for pid,(n,x) in sorted(c1.items()):
            S1V.append(n) 
            N1V.append(x) 
            CV.append(ci) 
        for pid,(n,x) in sorted(c2.items()):
            S2V.append(n) 
            N2V.append(x) 
    S1V = array(S1V)
    N1V = array(N1V)
    S2V = array(S2V)
    N2V = array(N2V)
    CV = array(CV)

    def vStats(S1V,N1V,S2V,N2V,CV):
        class VStats2:
                pass
        s = VStats2() 

        s.Ca = (CV==1).sum()
        s.Cu = (CV==0).sum()

        s.S1a = S1V[CV==1].sum()
        s.S1u = S1V[CV==0].sum() 
        s.N1a = N1V[CV==1].sum()
        s.N1u = N1V[CV==0].sum()
        s.SN1a = float(s.S1a) / s.N1a
        s.SN1u = float(s.S1u) / s.N1u

        s.S2a = S2V[CV==1].sum()
        s.S2u = S2V[CV==0].sum() 
        s.N2a = N2V[CV==1].sum()
        s.N2u = N2V[CV==0].sum()
        s.SN2a = float(s.S2a) / s.N2a
        s.SN2u = float(s.S2u) / s.N2u

        s.EN1a = s.SN1u * s.N1a
        s.EN2a = s.SN2u * s.N2a

        s.B = s.S1a + s.S2a
        s.EB = s.EN1a + s.EN2a
        s.delta = s.B - s.EB

        s.AD = 100. * s.delta / s.Ca 

        return s
    sReal = vStats(S1V,N1V,S2V,N2V,CV)
    r = [sReal]
    if nullI:
        r.append([vStats(S1V,N1V,S2V,N2V,permutation(CV)) for e in xrange(nullI)])
    if bootstrapI:
        bts = []
        for e in xrange(bootstrapI):
            ris = randint(len(CV),size=len(CV))
            bts.append(vStats(S1V[ris],N1V[ris],S2V[ris],N2V[ris],CV[ris]))
        r.append(bts)
    return r 


def empStats(sReal, sBckg, attribute):
    class EMPStats:
        pass 

    def one2twoPV(op):
        tp = 2.*op
        if tp > 1.0:
            tp = 2.*(1.-op)
        return tp

    es = EMPStats()
    es.v = getattr(sReal,attribute)
    es.rs = array([getattr(s,attribute) for s in sBckg])

    es.pvOne = (es.rs>=es.v).sum()/float(len(es.rs))
    es.pvTwo = one2twoPV(es.pvOne)

    rss = sort(es.rs)
    d = int(len(rss) * 0.025)
    es.left95 = rss[d]
    es.right95 = rss[-d]

    es.rsM = es.rs.mean()
    es.rsStd = es.rs.std()
    if es.rsStd == 0.0:
        es.z = nan
    else:
        es.z = (es.v - es.rsM) / es.rsStd
    es.pvOneAn = st.norm.sf(es.z)
    es.pvTwoAn = one2twoPV(es.pvOneAn)

    return es