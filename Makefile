.DELETE_ON_ERROR:

D=../src

all: Figure-1-cell-line-genetic-drift-filter.png \
    Table-1-small_scale.txt Table-2-allCNVs.txt Table-3-simplex-vs-multiplex.txt \
    Table-4-oneGeneCNVs.txt Table-5-smallscale-inter-coding_intronic.txt \
    SuppTable-1-smallscale-peripheral.txt \
    SuppFigure-1-CNVsByFilter.png SuppFigure-2-powerSnv.png \
    SuppFigure-3-parentalAges.png SuppFigure-4-powerCNV.png \
    SuppFigure-5-PCbyCNVGeneN.png SuppFigure-6-ratesVsAge.png \
    SuppFigure-7-ISB-power.png \
    SuppTable-2-propertyTable.txt

Figure-1-cell-line-genetic-drift-filter.png:
	python $D/figChildrenScatter.py $@

results.flag:
	python $D/result_tables.py . 1 && echo done > $@

Table-1-small_scale.txt: results.flag 
	python $D/mainTable-smallScale.py > $@ 

Table-2-allCNVs.txt: results.flag 
	python $D/mainTable-allCNVs.py > $@ 

Table-3-simplex-vs-multiplex.txt: results.flag 
	python $D/mainTable-simplex-multiplex-combined.py > $@ 

Table-4-oneGeneCNVs.txt: results.flag 
	python $D/mainTable-oneGeneCNVs.py > $@

Table-5-smallscale-inter-coding_intronic.txt: results.flag
	python $D/mainTable-intronic-peripheral.py inter-coding_intronic > $@

SuppTable-1-smallscale-peripheral.txt: results.flag
	python $D/mainTable-intronic-peripheral.py peripheral > $@

SuppFigure-1-CNVsByFilter.png: 
	python $D/figCNVsByFilter.py $@

SuppFigure-2-powerSnv.png:
	python $D/figSnvPower.py $@

SuppFigure-3-parentalAges.png: 
	python $D/figParentalAges.py $@

SuppFigure-4-powerCNV.png:
	python $D/figCNVPower.py $@

SuppFigure-5-PCbyCNVGeneN.png:
	python $D/figPCbyCNVgenes.py $@ 1

SuppFigure-6-ratesVsAge.png:
	python $D/figRatesVsAge.py $@

SuppFigure-7-ISB-power.png: results.flag
	python $D/figISBSignalPower.py $@ 1

SuppTable-2-propertyTable.txt:
	python $D/draw_float_property_hist.py . 2 8

