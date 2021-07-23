#!/usr/bin/env python

from diData import *
from collections import Counter, defaultdict
from scipy.stats import chi2_contingency
import scipy.stats as st

def count(chldrnIds,varFilterF):
    PSD = {pId:0 for pId in chldrnIds}
    for e in EVS:
        if not varFilterF(e): continue
        for pid in e.pids:
            if pid in PSD:
                PSD[pid] += 1
    return PSD 

def compareN2(cF1,cB1,cF2,cB2,nullI=0,bootstrapI=0,norStandard=None):
    assert cF1.keys() == cF2.keys()
    assert cB1.keys() == cB2.keys()

    N1V = []
    X1V = []
    N2V = []
    X2V = []
    SV = []
    for ci,(c1,c2) in enumerate([(cB1,cB2),(cF1,cF2)]):
        for pid,(n,x) in sorted(c1.items()):
            N1V.append(n) 
            X1V.append(x) 
            SV.append(ci) 
        for pid,(n,x) in sorted(c2.items()):
            N2V.append(n) 
            X2V.append(x) 
    N1V = array(N1V)
    X1V = array(X1V)
    N2V = array(N2V)
    X2V = array(X2V)
    SV = array(SV)

    def vStats(N1V,X1V,N2V,X2V,SV):
        class VStats:
                pass
        s = VStats() 

        s.A = (SV==1).sum()
        s.U = (SV==0).sum()

        s.N1a = N1V[SV==1].sum()
        s.N1u = N1V[SV==0].sum() 
        s.X1a = X1V[SV==1].sum()
        s.X1u = X1V[SV==0].sum()

        s.N2a = N2V[SV==1].sum()
        s.N2u = N2V[SV==0].sum() 
        s.X2a = X2V[SV==1].sum()
        s.X2u = X2V[SV==0].sum()

        s.R1a = float(s.N1a) / s.X1a
        s.R1u = float(s.N1u) / s.X1u

        s.R2a = float(s.N2a) / s.X2a
        s.R2u = float(s.N2u) / s.X2u

        s.EN1a = s.R1u * s.X1a
        s.EN2a = s.R2u * s.X2a

        s.B = s.N1a + s.N2a
        s.EB = s.EN1a + s.EN2a
        s.delta = s.B - s.EB

        s.AD = 100. * s.delta / s.A 

        return s
    sReal = vStats(N1V,X1V,N2V,X2V,SV)
    r = [sReal]
    if nullI:
        r.append([vStats(N1V,X1V,N2V,X2V,permutation(SV)) for e in xrange(nullI)])
    if bootstrapI:
        bts = []
        for e in xrange(bootstrapI):
            ris = randint(len(SV),size=len(SV))
            bts.append(vStats(N1V[ris],X1V[ris],N2V[ris],X2V[ris],SV[ris]))
        r.append(bts)
    return r 

def compareN(cF,cB,nullI=0,bootstrapI=0,norStandard=None):
    NV = []
    XV = []
    SV = []
    for ci,c in enumerate([cB,cF]):
        for pid,(n,x) in sorted(c.items()):
            NV.append(n) 
            XV.append(x) 
            SV.append(ci) 
    NV = array(NV)
    XV = array(XV)
    SV = array(SV)

    def vStats(NV,XV,SV):
        class VStats:
                pass
        s = VStats() 

        s.A = (SV==1).sum()
        s.U = (SV==0).sum()
        s.Na = NV[SV==1].sum()
        s.Nu = NV[SV==0].sum() 
        s.Xa = XV[SV==1].sum()
        s.Xu = XV[SV==0].sum()

        s.RTa = float(s.Na) / s.Xa
        s.RTu = float(s.Nu) / s.Xu

        s.ENa = s.RTu * s.Xa
        s.delta = s.Na - s.ENa

        s.AD = 100. * s.delta / s.A 

        if s.RTa:
            s.IR = 100. * (s.RTa-s.RTu) / s.RTa
        else:
            s.IR = 100. 

        if s.Na:
            s.CP = 100. * s.delta / s.Na 
        else:
            s.CP = 100. 

        s.RaR = float(s.Na)/s.A
        s.RuR = float(s.Nu)/s.U

        if norStandard:
            s.RaN = norStandard * s.Na/float(s.Xa)
            s.RuN = norStandard * s.Nu/float(s.Xu)
        return s
    sReal = vStats(NV,XV,SV)
    r = [sReal]
    if nullI:
        r.append([vStats(NV,XV,permutation(SV)) for e in xrange(nullI)])
    if bootstrapI:
        bts = []
        for e in xrange(bootstrapI):
            ris = randint(len(NV),size=len(NV))
            bts.append(vStats(NV[ris],XV[ris],SV[ris]))
        r.append(bts)
    return r 

def empStats(sReal, sBckg, a):
    class EMPStats:
        pass 

    def one2twoPV(op):
        tp = 2.*op
        if tp > 1.0:
            tp = 2.*(1.-op)
        return tp


    es = EMPStats()
    es.v = getattr(sReal,a)
    es.rs = array([getattr(s,a) for s in sBckg])

    es.pvOne = (es.rs>=es.v).sum()/float(len(es.rs))
    es.pvTwo = one2twoPV(es.pvOne)

    rss = sort(es.rs)
    d = int(len(rss) * 0.025)
    es.left95 = rss[d]
    es.right95 = rss[-d]

    es.rsM = es.rs.mean()
    es.rsStd = es.rs.std()
    es.z = (es.v - es.rsM) / es.rsStd
    es.pvOneAn = st.norm.sf(es.z)
    es.pvTwoAn = one2twoPV(es.pvOneAn)

    return es

def create_small_scale_result():
    CGG = [
        ["SSC unaffected", {pd.pId for pd in persons.values() if pd.coll == 'SSC' and pd.role == 'sib'}], 
        ["SSC affected", {pd.pId for pd in persons.values() if pd.coll == 'SSC' and pd.role == 'prb'}],
        ["AGRE affected", {pd.pId for pd in persons.values() if pd.coll == 'AGRE' and pd.role == 'prb'}]
    ]

    effTps = ['synonymous', 'LGD', 'nonsense', 'splice-site', 'frame-shift', 'missense']

    CNTS_N = {}

    SYN = count(persons, lambda e: (not e.location.startswith('chrX')) and e.eff == 'synonymous')
    for grpName,grpChIds in CGG:

        for effTp in effTps:

            if effTp == 'LGD':
                effTpS = set(['nonsense','frame-shift','splice-site'])
            else:
                effTpS = set([effTp])
            def supFltr(e):
                return (not e.location.startswith('chrX')) and e.eff in effTpS
            cnts = count(grpChIds, supFltr)
            CNTS_N[grpName,effTp] = {pid:(c,SYN[pid]) for pid,c in cnts.items()}

    def sm(g,e):
        return sum(CNTS[g,e].values())

    SSTF = open(outDir + '/small_scale_result_table.txt', 'w')
    ohcs = ['group','backgroundGroup', 'effect']
    ohcs += "sReal.U sReal.Nu sReal.Xu sReal.A sReal.Na sReal.Xa sReal.RaR sReal.ENa sReal.delta sReal.AD bcAD.left95 bcAD.right95 esAD.pvOne esAD.z esAD.pvOneAn sReal.IR bcIR.left95 bcIR.right95".split()

    SSTF.write("\t".join(ohcs)+'\n')
    globalRatios = {}
    obsNCntss = {}
    normSums = {} 
    for grpName,grpNameB in [
                ("SSC affected", "SSC unaffected"),
                ("AGRE affected","SSC unaffected"),
                ("SSC affected", "AGRE affected")]:
        for effTp in effTps:
            cntsB = CNTS_N[grpNameB,effTp]
            cntsF = CNTS_N[grpName,effTp]
            sReal, sNullBckg, sBtstrp = compareN(cntsF,cntsB,nullI=1000,bootstrapI=1000)
            esAD = empStats(sReal, sNullBckg, 'AD')
            bcAD = empStats(sReal, sBtstrp, 'AD')
            bcIR = empStats(sReal, sBtstrp, 'IR')

            cs = [grpName, grpNameB, effTp] 
    
            cs += map(str,[sReal.U, sReal.Nu, sReal.Xu, sReal.A, sReal.Na, sReal.Xa, sReal.RaR, sReal.ENa, sReal.delta, sReal.AD, bcAD.left95, bcAD.right95, esAD.pvOne, esAD.z, esAD.pvOneAn, sReal.IR, bcIR.left95, bcIR.right95])

            print "\t".join(map(str,cs))
            SSTF.write("\t".join(map(str,cs))+'\n')
    SSTF.close()

def create_merged_intron_result():
    def isIntercodingCNV(e):
        return e.vtype == 'CNV' and len(e.gns) and e.genomicRegion == 'inter-coding_intronic'


    indelVarTs= set(['del', 'ins'])
    def isTargetIntercodingIndel(e,genesP):
        return e.vtype in indelVarTs and e.genomicRegion == "inter-coding_intronic"  and e.gns & genesP 
    
    def isIntergenicIndel(e):
        return e.vtype in indelVarTs and e.eff == 'intergenic'

    def cntHelp(f):
        return count(persons, lambda e: (not e.location.startswith('chrX')) and f(e))

    ICNV = cntHelp(isIntercodingCNV)
    GIND = cntHelp(isIntergenicIndel)

    MIF = open(outDir + '/merged_intron_results.txt', 'w')
    ohcs = "genes geneNumber B EB delta AD AD.left AD.right AD.pval AD.z AD.pvOneAn".split()
    print "\t".join(ohcs)
    MIF.write("\t".join(ohcs) + "\n")

    targetDef = [
        'autism LGD',
        'all NDD LGD',
    ]
    for setK in targetDef:
        cn = 'number_of_' + CN(setK) + "_variants"
        genes = set(GENE['gene'][GENE[cn]>0])

        IIND = cntHelp(lambda e: isTargetIntercodingIndel(e,genes))

        CNTS_SEP_N = {}
        for role in ['prb','sib']:
            grpChIds = {pd.pId for pd in persons.values() if pd.coll == 'SSC' and pd.role == role}
            CNTS_SEP_N[role] = (
                        {pid:(ICNV[pid], 1.0) for pid in grpChIds},
                        {pid:(IIND[pid], GIND[pid]) for pid in grpChIds}
                    ) 

        cntsF_ICNV,cntsF_IINDs = CNTS_SEP_N['prb']
        cntsB_ICNV,cntsB_IINDs = CNTS_SEP_N['sib']
        sReal, sNullBckg, sBtstrp = compareN2(cntsF_ICNV, cntsB_ICNV, cntsF_IINDs, cntsB_IINDs, nullI=10000, bootstrapI=10000)
        esAD = empStats(sReal, sNullBckg, 'AD')
        bcAD = empStats(sReal, sBtstrp, 'AD')
        cs = map(str,[setK, len(genes), sReal.B, sReal.EB, sReal.delta, sReal.AD, bcAD.left95, bcAD.right95, esAD.pvOne, esAD.z, esAD.pvOneAn])
        print "\t".join(cs)
        MIF.write("\t".join(cs) + "\n")
    MIF.close()

def create_merged_smimplex_vs_multiplex_result():
    CGG = [
        ["SSC unaffected", {pd.pId for pd in persons.values() if pd.coll == 'SSC' and pd.role == 'sib'}], 
        ["SSC affected", {pd.pId for pd in persons.values() if pd.coll == 'SSC' and pd.role == 'prb'}],
        ["AGRE affected", {pd.pId for pd in persons.values() if pd.coll == 'AGRE' and pd.role == 'prb'}]
    ]
    def isCodingCNV(e):
        return e.vtype == 'CNV' and e.atts['size'] > 4000 and e.genomicRegion == 'coding'

    LGDEffs = set(['nonsense','frame-shift','splice-site'])
    def isLGD(e):
        return e.eff in LGDEffs
    
    def isSynonymous(e):
        return e.eff == 'synonymous'

    def cntHelp(f):
        return count(persons, lambda e: (not e.location.startswith('chrX')) and f(e))
    SYN = cntHelp(isSynonymous)
    LGD = cntHelp(isLGD)
    CNV = cntHelp(isCodingCNV)

    CNTS_SEP_N = {}
    for grpName,grpChIds in CGG:
        CNTS_SEP_N[grpName] = (
                    {pid:(CNV[pid], 1.0) for pid in grpChIds},
                    {pid:(LGD[pid],SYN[pid]) for pid in grpChIds}
                ) 

    MSMF = open(outDir + '/merged_simplex_multiplex_results.txt', 'w')

    ohcs = ["foreground group", "background group"]
    ohcs += "B EB delta AD AD.left AD.right AD.pval AD.z AD.pvOneAn".split()

    print "\t".join(ohcs)
    MSMF.write("\t".join(ohcs) + "\n")
    for grpF,grpB in [
                    ['SSC affected','SSC unaffected'],
                    ['AGRE affected','SSC unaffected'],
                    ['SSC affected','AGRE affected'] ]:

        cntsF_CNV,cntsF_LGDs = CNTS_SEP_N[grpF]
        cntsB_CNV,cntsB_LGDs = CNTS_SEP_N[grpB]
        sReal, sNullBckg, sBtstrp = compareN2(cntsF_CNV, cntsB_CNV, cntsF_LGDs, cntsB_LGDs ,nullI=1000,bootstrapI=1000)
        esAD = empStats(sReal, sNullBckg, 'AD')
        bcAD = empStats(sReal, sBtstrp, 'AD')
        cs = [grpF, grpB]
        cs += map(str,[sReal.B, sReal.EB, sReal.delta, sReal.AD, bcAD.left95, bcAD.right95, esAD.pvOne, esAD.z, esAD.pvOneAn])
        print "\t".join(map(str,cs))
        MSMF.write("\t".join(map(str,cs)) + "\n")
    MSMF.close()

def create_CNV_result():
    CGG = [
        ["SSC unaffected", {pd.pId for pd in persons.values() if pd.coll == 'SSC' and pd.role == 'sib'}], 
        ["SSC affected", {pd.pId for pd in persons.values() if pd.coll == 'SSC' and pd.role == 'prb'}],
        ["AGRE affected", {pd.pId for pd in persons.values() if pd.coll == 'AGRE' and pd.role == 'prb'}]
    ]

    CRTF = open(outDir + '/CNV_result_table.txt', 'w')

    for grpName,grpChIds in CGG:
        print grpName, len(grpChIds)


    def okCNV(e):
        return e.atts['size'] > 4000 and e.vtype == 'CNV'
        

    allEffTps = [
        ['all', lambda e: okCNV(e)],
        ['coding', lambda e: okCNV(e) and e.genomicRegion == 'coding'],
        ['intergenic', lambda e: okCNV(e) and e.genomicRegion == 'intergenic'],
        ['genic noncoding', lambda e: okCNV(e) and e.genomicRegion != 'coding' and e.genomicRegion != 'intergenic']
    ]

    def oneGeneCNV(e):
        return e.vtype == 'CNV' and len(e.gns) == 1


    oneGeneEffTps = [
        ['all',                  lambda e: oneGeneCNV(e) ],
        ['coding',               lambda e: oneGeneCNV(e) and e.genomicRegion == 'coding'],
        ['intercoding intronic', lambda e: oneGeneCNV(e) and e.genomicRegion == 'inter-coding_intronic'],
        ['peripheral',           lambda e: oneGeneCNV(e) and e.genomicRegion == 'peripheral']
    ]

    allVTps = [
        ['deletion',             lambda e: e.vtype == "CNV" and e.eff == 'deletion'],
        ['duplication',          lambda e: e.vtype == "CNV" and e.eff == 'duplication'] ,
        ['cnvs',                 lambda e: e.vtype == "CNV" ]
    ]

    sectionDef =  zip(['ALL EVENTS','ONE CODING GENE'],[allEffTps,oneGeneEffTps])

    CNTS = {}
    CNVs = [e for e in EVS if e.vtype=="CNV"]
    def countCNVs(chldrnIds,varFilterF):
        PSD = {pId:0 for pId in chldrnIds}
        for e in CNVs:
            if not varFilterF(e): continue
            for pid in e.pids:
                if pid in PSD:
                    PSD[pid] += 1
        return {pid:(c,1) for pid,c in PSD.items()}
        
    for grpName,grpChIds in CGG:
        for secName,effTps in sectionDef:
            for effTp,effFF in effTps:
                for vTp,vtFF in allVTps:
                    CNTS[grpName,secName,effTp,vTp] = countCNVs(grpChIds, lambda e: effFF(e) and vtFF(e))

    ohcs = ["section","effect", "variant", "group", "backgroundGroup"]
    ohcs += "sReal.U sReal.Nu sReal.Xu sReal.A sReal.Na sReal.Xa sReal.RaR sReal.ENa sReal.delta sReal.AD bcAD.left95 bcAD.right95 esAD.pvOne esAD.z esAD.pvOneAn sReal.IR bcIR.left95 bcIR.right95".split()
    CRTF.write("\t".join(ohcs) + "\n")
    for secName,effTps in sectionDef:
        for effTp,effFF in effTps:
            for vTp,vtFF in allVTps:
                for grpName,grpNameB in [
                            ("SSC affected", "SSC unaffected"),
                            ("AGRE affected","SSC unaffected"),
                            ("SSC affected", "AGRE affected")]:
                    cntsF = CNTS[grpName, secName,effTp,vTp]
                    cntsT = CNTS[grpNameB,secName,effTp,vTp]

                    sReal, sNullBckg, sBtstrp = compareN(cntsF,cntsT,nullI=1000,bootstrapI=1000)
                    esAD = empStats(sReal, sNullBckg, 'AD')
                    bcAD = empStats(sReal, sBtstrp, 'AD')
                    bcIR = empStats(sReal, sBtstrp, 'IR')
                    cs = [secName, effTp, vTp, grpName, grpNameB] 
                    cs += map(str,[sReal.U, sReal.Nu, sReal.Xu, sReal.A, sReal.Na, sReal.Xa, sReal.RaR, sReal.ENa, sReal.delta, sReal.AD, bcAD.left95, bcAD.right95, esAD.pvOne, esAD.z, esAD.pvOneAn, sReal.IR, bcIR.left95, bcIR.right95])
                    CRTF.write("\t".join(cs) + "\n")

    for grpName,grpChIds in CGG:
        CRTF.write("#\t%s: %d\n" % (grpName, len(grpChIds)))
    CRTF.close()

def create_intronic_result():

    prbs = {pd.pId for pd in persons.values() if pd.coll == 'SSC' and pd.role == 'prb'}
    sibs = {pd.pId for pd in persons.values() if pd.coll == 'SSC' and pd.role == 'sib'}

    mtSets = [
        'all genes',
        'autism LGD',
        'autism missense',
        'autism synonymous',
        'all NDD LGD',
        'all NDD missense',
        'all NDD synonymous',
    ]

    IRTF = open(outDir + '/intronic_result_table.txt', 'w')

    hcs = "set eventType intronType setGeneNumber".split(" ")
    hcs += "sReal.U sReal.Nu sReal.Xu sReal.A sReal.Na sReal.Xa sReal.RaR sReal.ENa sReal.delta sReal.AD bcAD.left95 bcAD.right95 esAD.pvOne esAD.z esAD.pvOneAn sReal.IR bcIR.left95 bcIR.right95".split()
    IRTF.write("\t".join(hcs)+"\n")
    print "\t".join(hcs)

    normVariant = {
        "sub":   count(persons, lambda e: (not e.location.startswith('chrX')) and e.genomicRegion == 'intergenic' and e.vtype in ['sub'] and len(e.pids) == 1),
        "indel": count(persons, lambda e: (not e.location.startswith('chrX')) and e.genomicRegion == 'intergenic' and e.vtype in ['ins','del'] and len(e.pids) == 1)
    }

    for stK in mtSets:
        if stK == 'all genes':
            iis = ones(len(GENE),dtype=bool)
        else:
            cn = 'number_of_' + CN(stK) + "_variants"
            iis = GENE[cn] > 0
        
        for eT,varTs in zip(["indel","sub"],[set(["del","ins"]),set(["sub"])]):

            for effT in ["inter-coding_intronic","peripheral"]:
                
                gnsS = set(GENE['gene'][iis])

                cs = [stK,eT,effT,iis.sum()] 
                def eventFF(e):
                    return (not e.location.startswith('chrX')) and e.vtype in varTs and e.genomicRegion==effT and (e.gns & gnsS) and e.eff != "splice-site" and len(e.pids) == 1
                def countHelper(chldrn):
                    cnt = count(chldrn,eventFF)
                    return {pid:(v,normVariant[eT][pid]) for pid,v in cnt.items()}
                prbCs = countHelper(prbs)
                sibCs = countHelper(sibs)

                sReal, sNullBckg, sBtstrp = compareN(prbCs,sibCs,nullI=1000,bootstrapI=1000)
                esAD = empStats(sReal, sNullBckg, 'AD')
                bcAD = empStats(sReal, sBtstrp, 'AD')
                bcIR = empStats(sReal, sBtstrp, 'IR')
                cs += map(str,[sReal.U, sReal.Nu, sReal.Xu, sReal.A, sReal.Na, sReal.Xa, sReal.RaR, sReal.ENa, sReal.delta, sReal.AD, bcAD.left95, bcAD.right95, esAD.pvOne, esAD.z, esAD.pvOneAn, sReal.IR, bcIR.left95, bcIR.right95])
                print "\t".join(map(str,cs))
                IRTF.write("\t".join(map(str,cs))+"\n")

    IRTF.close()

if __name__ == "__main__":
    tb = 'all' if len(sys.argv) < 2 else sys.argv[1]
    outDir = '.' if len(sys.argv) < 3 else sys.argv[2]
    seedV = None if len(sys.argv) < 4 else int(sys.argv[3])

    if seedV: seed(seedV)
        
    EVS = loadEVS(['De novo CNV in SSC and AGRE', 'Small scale de novo in SSC', 'Small scale de novo in AGRE'])

    if tb == 'all':
        create_small_scale_result()
        create_intronic_result()
        create_CNV_result()
        create_merged_smimplex_vs_multiplex_result()
        create_merged_intron_result()
    elif tb == 'small_scale':
        create_small_scale_result()
    elif tb == 'CNV':
        create_CNV_result()
    elif tb == 'intronic':
        create_intronic_result()
    elif tb == 'merged_simplex_multiplex':
        create_merged_smimplex_vs_multiplex_result()
    elif tb == 'merged_intron':
        create_merged_intron_result()
    else:
        print "The argument should be all or small_scale or CNV or intronic or merged_simplex_multiplex or merged_intron" 
    
