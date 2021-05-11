D=./src

all: Figure-1-genomic-drift-outliers.png \
    Table-1-small_scale.txt Table-2-allCNVs.txt Table-3-simplex-vs-multiplex.txt \
    Table-4-oneGeneCNVs.txt Table-5-smallscale-inter-coding_intronic.txt \
    SuppTable-6-smallscale-peripheral.txt \
    SuppFigure-1-CNVsByFilter.png SuppFigure-2-powerSnv.png \
    SuppFigure-3-parentalAges.png SuppFigure-4-powerCNV.png \
    SuppFigure-5-PCbyCNVGeneN.png SuppFigure-6-ratesVsAge.png \
    SuppFigure-7-ISB-power.png
 
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
	$D/figPCbyCNVgenes.py $@

SuppFigure-6-ratesVsAge.png:
	$D/figRatesVsAge.py $@

SuppFigure-7-ISB-power.png:
	$D/figISBSignalPower.py $@


propPNGs.flag: $D/draw_float_property_hist.py 
	(time (mkdir -p propPNGs &&  $< propPNGs  > log/$@-out.txt 2> log/$@-err.txt && touch $@)) 2> log/$@-time.txt


ss-model-figures.flag: $D/splice_site_model.py ss-model-data.flag
	(time $<    > log/$@-out.txt 2> log/$@-err.txt && touch $@) 2> log/$@-time.txt

table_sanders_missing.txt: $D/table_sanders_missing.py ${DATA_FILES}
	(time $<    > log/$@-out.txt 2> log/$@-err.txt && touch $@) 2> log/$@-time.txt


