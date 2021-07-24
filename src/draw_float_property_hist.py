#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')

import sys
from pylab import *
from scipy.stats import ttest_ind,ks_2samp,ranksums,chi2_contingency,fisher_exact,mannwhitneyu
from scipy.stats import rankdata
from diData import loadEVS, GENE, CN, persons
from pV2Str import pV2Str

outDir = "."
tableNumber = 1
startFigureNumber = 1

if len(sys.argv)>1:
    outDir = sys.argv[1]
if len(sys.argv)>2:
    tableNumber = int(sys.argv[2])
if len(sys.argv)>3:
    startFigureNumber = int(sys.argv[3])


gnsS = set(GENE['gene'][GENE[CN('number of all NDD LGD variants')]>0])
EVS = loadEVS(['Small scale de novo in SSC'])

propDefs = [
#    prop,                                   doAbs, legendLeft, bestRank 
    ("variant size",                         True , False,      'max'    ),
    ("intron length",                        False, False,      'min'    ),
    ("distance from splice-site",            True , False,      'min'    ),
    ("ORF length",                           False, False,      'max'    ),

    ('__SECTION__',         'SpliceAI Scores',                    '', '' ),
    ("SpliceAI DS_AG score",                 False, False,      'max'    ),
    ("SpliceAI DS_AL score",                 False, False,      'max'    ),
    ("SpliceAI DS_DG score",                 False, False,      'max'    ),
    ("SpliceAI DS_DL score",                 False, False,      'max'    ),
    ("SpliceAI MAX_DS score",                False, False,      'max'    ),

    ('__SECTION__',         'Simple Splice Model Scores',         '', '' ),
    ("acceptor 'alt' score",                 False, False,      'max'    ),
    ("acceptor 'ref' score",                 False, False,      'max'    ),
    ("acceptor 'alt-ref' score",             False, False,      'max'    ),
    ("donor 'alt' score",                    False, False,      'max'    ),
    ("donor 'ref' score",                    False, True,       'max'    ),
    ("donor 'alt-ref' score",                False, False,      'max'    ),

    ('__SECTION__',         'Conservation Scores',                '', '' ),
    ("phylop, 100 vertebrates score",        False, False,      'max'    ),
    ("phylop, 30 vertebrates score",         False, True,       'max'    ),
    ("phylop, 20 vertebrates score",         False, True,       'max'    ),
    ("phylop, 7 vertebrates score",          False, True,       'max'    ),
    ("phastCons, 100 vertebrates score",     False, False,      'max'    ),
    ("phastCons, 30 vertebrates score",      False, False,      'max'    ),
    ("phastCons, 20 vertebrates score",      False, False,      'max'    ),
    ("phastCons, 7 vertebrates score",       False, False,      'max'    ),
    
    ("CADD score",                           False, False,      'max'    ),

    ('__SECTION__',         'Min Rank Aggregated score',          '', '' ),
    ("minimum property rank",                False, False,      'min'    ),

]

clps = [clp for varT in ["indel","sub"] for clp in [((varT,"T","P"),(varT,"T","S")), ((varT,"A","P"),(varT,"A","S"))]]

def clS(cl):
    return "%s (%s,%s)" % ( \
        "IID" if cl[0]=='indel' else 'ISB', \
        "target genes" if cl[1]=="T" else "all genes", \
        "affected" if cl[2]=="P" else "unaffected")

    return "".join([x[0] for x in cl])

def getEvents(vrTps,skipShared=True,quadsOnly=True,skipX=True):
    for e in EVS:
        if e.coll != "SSC"                                : continue
        if e.vtype not in vrTps                           : continue
        if e.eff == 'splice-site'                         : continue
        if e.genomicRegion != 'inter-coding_intronic'     : continue
        if skipX and e.location.startswith('chrX')        : continue
        if skipShared and len(e.pids)>1                   : continue
        yield e

def grankByProp(evs,propN,doAbs,bestRank):
    assert bestRank in ['min', 'max']
    vs = array([e.atts[CN(propN)] for e in evs])

    msngInds = isnan(vs)

    if doAbs: 
        vs = abs(vs)

    kvs = vs[logical_not(msngInds)]
    if bestRank == 'min':
        vs[msngInds] = max(kvs) + 1
    else:
        vs[msngInds] = min(kvs) - 1

    if bestRank == 'max':
        vs = -vs
    rnks = rankdata(vs)
    return rnks

def minRanks(evs):
    allRnks = vstack([ grankByProp(evs,propN,doAbs,bestRank)  
            for propN,doAbs,legendLeft,bestRank in propDefs 
                if propN != 'minimum property rank' and propN != '__SECTION__']).T
    return allRnks.min(axis=1)

TF = open(outDir + "/SuppTable-%s-propertyTable.txt" % (tableNumber),"w")
TF.write("\t".join(['property','Supp. Fig. N.'] + ["-".join(map(clS,clp)) for clp in clps]) + "\n")

propI = 0
for propN,doAbs,legendLeft,bestRank in propDefs:
    if propN == '__SECTION__':
        TF.write(doAbs + "\n")
        continue

    propCN = CN(propN) 

    vsByClass = {}
    mns = []
    mxs = []

    for varT in ["indel","sub"]:
        vrTps = set(["del","ins"]) if varT == 'indel' else set([varT])
        values_empty=False
        
        vEvs = array(list(getEvents(vrTps)))

        if propN == "minimum property rank":
            prop = minRanks(vEvs)
        else:
            prop = array([e.atts[propCN] for e in vEvs])

        vEvs = vEvs[isfinite(prop)] 
        prop = prop[isfinite(prop)] 

        if doAbs:
            prop = abs(prop)

        try:
            mns.append(min(prop)) 
            mxs.append(max(prop)) 
        except ValueError:
            values_empty=True  
            print("no values for ", propN)

        inPrb = array(['prb' in {persons[pid].role for pid in e.pids} for e in vEvs])
        inSib = array(['sib' in {persons[pid].role for pid in e.pids} for e in vEvs])

        inSet = array([bool(set(e.gns & gnsS)) for e in vEvs],dtype=bool)

        inPrbS = logical_and(inPrb,inSet)
        inSibS = logical_and(inSib,inSet)
        
        if len(prop) != 0:
            vsByClass[(varT,'A','P')] = prop[inPrb]
            vsByClass[(varT,'A','S')] = prop[inSib]
            vsByClass[(varT,'T','P')] = prop[inPrbS]
            vsByClass[(varT,'T','S')] = prop[inSibS]   
      
    if len(mns) != 0:
        mn = min(mns)
    if len(mxs) != 0:
        mx = max(mxs)

    supFN = startFigureNumber + propI
    propI += 1

    clf()
    tbCs = [propN,str(supFN)]
    for clpi,(cl1,cl2) in enumerate(clps):
        bns = 30
        subplot(2,2,clpi+1)

        if len(vsByClass) != 0:
            vs1 = vsByClass[cl1]
            vs2 = vsByClass[cl2]
        else:
            vs1=[]
            vs2=[]
        a,bns,b = hist(vs1,bins=bns,alpha=0.5,density=True,color='blue')
        a,bns,b = hist(vs2,bins=bns,alpha=0.5,density=True,color='green')

        if len(vs1) != 0 and len(vs2) != 0:
            aKs = ks_2samp(vs1,vs2)[1] 
            aRs = ranksums(vs1,vs2)[1] 
            try:
                aMs = mannwhitneyu(vs1,vs2,alternative='two-sided')[1]
            except ValueError:
                aMs = np.nan
            aTt = ttest_ind(vs1,vs2)[1]

            pref = ''
            if clpi in [0,1]:
                pref = "\n\n"
            tto = title(pref + "pVals: mannW: %s, ks: %s, ttest: %s" % (pV2Str(aMs),pV2Str(aKs),pV2Str(aTt)))
            if clpi == 0:
                gcf().text(0.3,0.97,'target genes',fontsize=16,ha='center')
            if clpi == 1:
                gcf().text(0.76,0.97,'all genes',fontsize=16,ha='center')

            if isnan(aMs):
                tbCs.append('')
            else:
                tbCs.append(pV2Str(aMs))


            if clpi == 0:
                ylabel('de novo inter-coding intronic\nindels (IID)',fontsize=16)
            elif clpi == 2:
                ylabel('de novo inter-coding intronic\nsubstitutions (ISB)',fontsize=16)
            gca().get_yaxis().set_ticks([])
            if doAbs:
                xlabel("abs(" + propN + ")")
            else:
                xlabel(propN)

            def lgndS(cl,vs):
                return "%s (n: %d, m: %.3f, s: %.2f )" % \
                        ("affected    " if cl[2]=="P" else "unaffected", \
                         len(vs),vs.mean(), vs.std())
     
            lgnd =  [lgndS(cl1,vs1),lgndS(cl2,vs2)]

            if legendLeft:
                legend(lgnd,loc="upper left")
            else:
                legend(lgnd)
    
    gcf().set_size_inches(16,10)
    tight_layout()
    gcf().savefig("%s/SuppFigure-%d-hists-%s.png" % (outDir,supFN,propCN))

    TF.write("\t".join(tbCs) + "\n") 
TF.close()
