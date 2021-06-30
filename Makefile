.DELETE_ON_ERROR:

D=../src

all: Figure-1-genomic-drift-outliers.png \
    Table-1-small_scale.txt Table-2-allCNVs.txt Table-3-simplex-vs-multiplex.txt \
    Table-4-oneGeneCNVs.txt Table-5-smallscale-inter-coding_intronic.txt \
    SuppTable-6-smallscale-peripheral.txt \
    SuppFigure-1-CNVsByFilter.png SuppFigure-2-powerSnv.png \
    SuppFigure-3-parentalAges.png SuppFigure-4-powerCNV.png \
    SuppFigure-5-PCbyCNVGeneN.png SuppFigure-6-ratesVsAge.png \
    SuppFigure-7-ISB-power.png \
    SuppTable-8-propertyTable.txt

Figure-1-genomic-drift-outliers.png:
	$D/figChildrenScatter.py $@

results.flag: 
	$D/result_tables.py all . 1 && touch $@

Table-1-small_scale.txt: results.flag 
	$D/mainTable-smallScale.py > $@ 

Table-2-allCNVs.txt: results.flag 
	$D/mainTable-allCNVs.py > $@ 

Table-3-simplex-vs-multiplex.txt: results.flag 
	$D/mainTable-simplex-multiplex-combined.py > $@ 

Table-4-oneGeneCNVs.txt: results.flag 
	$D/mainTable-oneGeneCNVs.py > $@

Table-5-smallscale-inter-coding_intronic.txt: results.flag
	$D/mainTable.py inter-coding_intronic > $@

SuppTable-6-smallscale-peripheral.txt: results.flag
	$D/mainTable.py peripheral > $@

SuppFigure-1-CNVsByFilter.png: 
	$D/figCNVsByFilter.py $@

SuppFigure-2-powerSnv.png:
	$D/figSnvPower.py $@

SuppFigure-3-parentalAges.png: 
	$D/figParentalAges.py $@

SuppFigure-4-powerCNV.png:
	$D/figCNVPower.py $@

SuppFigure-5-PCbyCNVGeneN.png:
	$D/figPCbyCNVgenes.py $@ 1

SuppFigure-6-ratesVsAge.png:
	$D/figRatesVsAge.py $@

SuppFigure-7-ISB-power.png: results.flag
	$D/figISBSignalPower.py $@ 1

SuppTable-8-propertyTable.txt:
	$D/draw_float_property_hist.py

