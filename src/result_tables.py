#!/usr/bin/env python

from pylab import *
from collections import Counter, defaultdict
from scipy.stats import chi2_contingency
import scipy.stats as st
from diData import GENE, persons, CN, loadEVS
from methods import compare_subject_variant_class_in_two_groups_using_normalization_variant_class, empStats
from methods import compare_joinly_two_subject_variant_classes_in_two_groups_using_the_separate_normalization_variant_classes

_EVS = loadEVS(['De novo CNV in SSC and AGRE', 'Small scale de novo in SSC', 'Small scale de novo in AGRE'])

def count(chldrnIds,varFilterF):
    PSD = {pId:0 for pId in chldrnIds}
    for e in _EVS:
        if not varFilterF(e): continue
        for pid in e.pids:
            if pid in PSD:
                PSD[pid] += 1
    return PSD 

_CGG = [
    ["SSC unaffected", {pd.pId for pd in persons.values() if pd.coll == 'SSC' and pd.role == 'sib'}], 
    ["SSC affected", {pd.pId for pd in persons.values() if pd.coll == 'SSC' and pd.role == 'prb'}],
    ["AGRE affected", {pd.pId for pd in persons.values() if pd.coll == 'AGRE' and pd.role == 'prb'}]
]

_oneSubjectClassStatsColumns = ("sReal.Cu sReal.Su sReal.Nu " + \
                               "sReal.Ca sReal.Sa sReal.Na " + \
                               "sReal.RSa sReal.ESa sReal.delta " + \
                               "sReal.AD bcAD.left95 bcAD.right95 esAD.pvOne esAD.z esAD.pvOneAn " + \
                               "sReal.PC bcPC.left95 bcPC.right95").split()

def _computeOneSubjectStatsColumns(countsInAffected,countsInUnaffected):
    sReal, sNullBckg, sBtstrp = compare_subject_variant_class_in_two_groups_using_normalization_variant_class(countsInAffected,countsInUnaffected,nullI=1000,bootstrapI=1000)
    esAD = empStats(sReal, sNullBckg, 'AD')
    bcAD = empStats(sReal, sBtstrp, 'AD')
    bcPC = empStats(sReal, sBtstrp, 'PC')  
    return map(str,[sReal.Cu, sReal.Su, sReal.Nu, sReal.Ca, sReal.Sa, sReal.Na, sReal.RSa, sReal.ESa, sReal.delta, sReal.AD, bcAD.left95, bcAD.right95, esAD.pvOne, esAD.z, esAD.pvOneAn, sReal.PC, bcPC.left95, bcPC.right95])

def create_small_scale_result():
    effTps = ['synonymous', 'LGD', 'nonsense', 'splice-site', 'frame-shift', 'missense']

    CNTS_N = {}

    SYN = count(persons, lambda e: (not e.location.startswith('chrX')) and e.eff == 'synonymous')
    for grpName,grpChIds in _CGG:
        for effTp in effTps:
            if effTp == 'LGD':
                effTpS = set(['nonsense','frame-shift','splice-site'])
            else:
                effTpS = set([effTp])
            def supFltr(e):
                return (not e.location.startswith('chrX')) and e.eff in effTpS
            cnts = count(grpChIds, supFltr)
            CNTS_N[grpName,effTp] = {pid:(c,SYN[pid]) for pid,c in cnts.items()}

    SSTF = open(outDir + '/small_scale_result_table.txt', 'w')
    ohcs = ['group','backgroundGroup', 'effect']
    ohcs += _oneSubjectClassStatsColumns

    SSTF.write("\t".join(ohcs)+'\n')
    for grpName,grpNameB in [
                ("SSC affected", "SSC unaffected"),
                ("AGRE affected","SSC unaffected"),
                ("SSC affected", "AGRE affected")]:
        for effTp in effTps:
            cntsF = CNTS_N[grpName,effTp]
            cntsB = CNTS_N[grpNameB,effTp]

            cs = [grpName, grpNameB, effTp] + \
                 _computeOneSubjectStatsColumns(cntsF,cntsB)

            SSTF.write("\t".join(map(str,cs))+'\n')
    SSTF.close()

def create_CNV_result():
    CRTF = open(outDir + '/CNV_result_table.txt', 'w')

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
    CNVs = [e for e in _EVS if e.vtype=="CNV"]
    def countCNVs(chldrnIds,varFilterF):
        PSD = {pId:0 for pId in chldrnIds}
        for e in CNVs:
            if not varFilterF(e): continue
            for pid in e.pids:
                if pid in PSD:
                    PSD[pid] += 1
        return {pid:(c,1) for pid,c in PSD.items()}
        
    for grpName,grpChIds in _CGG:
        for secName,effTps in sectionDef:
            for effTp,effFF in effTps:
                for vTp,vtFF in allVTps:
                    CNTS[grpName,secName,effTp,vTp] = countCNVs(grpChIds, lambda e: effFF(e) and vtFF(e))

    ohcs = ["section","effect", "variant", "group", "backgroundGroup"]
    ohcs += _oneSubjectClassStatsColumns
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

                    cs = [secName, effTp, vTp, grpName, grpNameB] + \
                         _computeOneSubjectStatsColumns(cntsF,cntsT)

                    CRTF.write("\t".join(cs) + "\n")

    for grpName,grpChIds in _CGG:
        CRTF.write("#\t%s: %d\n" % (grpName, len(grpChIds)))
    CRTF.close()

def create_merged_smimplex_vs_multiplex_result():
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
    for grpName,grpChIds in _CGG:
        CNTS_SEP_N[grpName] = (
                    {pid:(CNV[pid], 1.0) for pid in grpChIds},
                    {pid:(LGD[pid],SYN[pid]) for pid in grpChIds}
                ) 

    MSMF = open(outDir + '/merged_simplex_multiplex_results.txt', 'w')

    ohcs = ["foreground group", "background group"]
    ohcs += "B EB delta AD AD.left AD.right AD.pval AD.z AD.pvOneAn".split()

    MSMF.write("\t".join(ohcs) + "\n")
    for grpF,grpB in [
                    ['SSC affected','SSC unaffected'],
                    ['AGRE affected','SSC unaffected'],
                    ['SSC affected','AGRE affected'] ]:

        cntsF_CNV,cntsF_LGDs = CNTS_SEP_N[grpF]
        cntsB_CNV,cntsB_LGDs = CNTS_SEP_N[grpB]
        sReal, sNullBckg, sBtstrp = compare_joinly_two_subject_variant_classes_in_two_groups_using_the_separate_normalization_variant_classes(cntsF_CNV, cntsB_CNV, cntsF_LGDs, cntsB_LGDs ,nullI=1000,bootstrapI=1000)
        esAD = empStats(sReal, sNullBckg, 'AD')
        bcAD = empStats(sReal, sBtstrp, 'AD')
        cs = [grpF, grpB]
        cs += map(str,[sReal.B, sReal.EB, sReal.delta, sReal.AD, bcAD.left95, bcAD.right95, esAD.pvOne, esAD.z, esAD.pvOneAn])
        MSMF.write("\t".join(map(str,cs)) + "\n")
    MSMF.close()

def create_intronic_result():
    CGGD = {k:iis for k,iis in _CGG}
    prbs = CGGD["SSC affected"]
    sibs = CGGD["SSC unaffected"]
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

    hcs = "geneSet eventType regionType geneSetGeneNumber".split(" ")
    hcs += _oneSubjectClassStatsColumns
    IRTF.write("\t".join(hcs)+"\n")

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

                def eventFF(e):
                    return (not e.location.startswith('chrX')) and e.vtype in varTs and e.genomicRegion==effT and (e.gns & gnsS) and e.eff != "splice-site" and len(e.pids) == 1
                def countHelper(chldrn):
                    cnt = count(chldrn,eventFF)
                    return {pid:(v,normVariant[eT][pid]) for pid,v in cnt.items()}

                prbCs = countHelper(prbs)
                sibCs = countHelper(sibs)

                cs = [stK,eT,effT,iis.sum()] + \
                    _computeOneSubjectStatsColumns(prbCs,sibCs)

                IRTF.write("\t".join(map(str,cs))+"\n")
    IRTF.close()

if __name__ == "__main__":
    outDir = '.' if len(sys.argv) < 2 else sys.argv[1]
    seedV = None if len(sys.argv) < 3 else int(sys.argv[2])

    if seedV: seed(seedV)
    
    create_small_scale_result()
    create_intronic_result()
    create_CNV_result()
    create_merged_smimplex_vs_multiplex_result()