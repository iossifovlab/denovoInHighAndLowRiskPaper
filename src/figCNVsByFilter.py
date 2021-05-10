#!/usr/bin/env python

from diData import * 
from collections import Counter

CNVR = genfromtxt(resDir + "/CNV_raw_table.txt", delimiter='\t',dtype=None,names=True, case_sensitive=True,encoding=None)

okChldrnN = {pd.pId:0 for pd in persons.values()}
flChldrnN = {pd.pId:0 for pd in outlierPersons.values()}

assert not set(okChldrnN) & set(flChldrnN)

ratios = defaultdict(list) 
for cnv in CNVR:
    # filter out the unique high quality variants
    if not cnv['if_merged']: continue
    if cnv['parents_dels'] == 'NA' or cnv['parents_dups'] == 'NA' or int(cnv['parents_dels']) > 5 or int(cnv['parents_dups']) > 5: continue
    if cnv['ploidy_measure'] == 'NA' or float(cnv['ploidy_measure']) < 0.85 or float(cnv['ploidy_measure']) > 1.15: continue
    # print 'ha'

    pid = cnv['sampleID']
    if pid in okChldrnN:
        okChldrnN[pid] += 1
        assert cnv['in_summary']
    elif pid in flChldrnN:
        flChldrnN[pid] += 1
        assert not cnv['in_summary']
    else:
        print "BREH: unknown pid", pid
        assert not cnv['in_summary']
        continue

    ratioS = cnv['region_median_ratio_for_the_child']
    if ratioS != 'NA': 
        ratio = float(ratioS)
        ok = pid in okChldrnN
        ratios[ok,cnv['polarity']].append(ratio)
clf()
xs = range(1,max(max(okChldrnN.values()),max(flChldrnN.values()))+1)
ttls = ['%d children that pass the SNV filter' % len(okChldrnN),
        '%d children that fail the SNV filter' % len(flChldrnN)]
spn = 1
for ni,N in enumerate([okChldrnN,flChldrnN]):
    subplot(2,2,spn)
    spn += 1
    cs = Counter(N.values())
    ps = [100.*float(cs[x])/len(N) for x in xs] 
    # print ni,cs,ns
    bar(xs,ps)
    ylim([0,15])
    ylabel('percent of children with\nx de novo CNVs')
    xlabel('number of CNVs')
    title(ttls[ni])

for plrty in ['deletion','duplication']:
    subplot(2,2,spn)
    spn += 1
    bns = 10 
    a,bns,c = hist(ratios[True,plrty], bins=bns,density=True,alpha=0.5,label='ok (%d %ss)'       % (len(ratios[True,plrty]), plrty))
    a,bns,c = hist(ratios[False,plrty],bins=bns,density=True,alpha=0.5,label='filtered (%d %ss)' % (len(ratios[False,plrty]),plrty))    
    legend()
    xlabel('CNV ratio')
    ylabel('pdf')
    title('CNV ratios for ' + plrty)

gcf().set_size_inches(10,8)
tight_layout()
gcf().savefig("CNVsByFilter.png")

# show()
