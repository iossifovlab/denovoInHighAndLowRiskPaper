SHELL=/bin/bash -o pipefail
.DELETE_ON_ERROR:

D=./src
I=./input


all: results.flag \
    main_table_nonCoding_peripheral.txt main_table_nonCoding_inter-coding_intronic.txt \
    main_table_small_scale.txt mainTable-simplex-multiplex-combined.txt \
    main_table_allCNVs.txt main_table_oneGeneCNVs.txt

#    children_scatter.png genomic_region_terms.png results.flag  \
#    CDintron-list.txt ss-model-data.flag ss-model-figures.flag \
#    SSC-children-trioDenovo-stats.txt AGRE-children-trioDenovo-stats.txt \
#    powerCNV.png powerSnv.png CNVsByFilter.png figParentalAges.png IRbyCNVGeneN.png \
#    figRatesVsAge.png ISB-power.png \
#    SSC_small_ann_funProp_denovo_table.txt propPNGs.flag \
#    AGRE_small_ann_funProp_denovo_table.txt suppTables.flag \
#    hmm_figure.png table_sanders_missing.txt

# clean:
# 	rm -f children_table.txt ${DN_TBLS_RAW} ${DN_TBLS_ANN} ssc_intronic_res_table.txt

results.flag: $D/result_tables.py ${DATA_FILES}
	(time $< all > log/$@-out.txt 2> log/$@-err.txt && touch $@) 2> log/$@-time.txt

main_table_small_scale.txt: results.flag $D/mainTable-smallScale.py
	(time $D/mainTable-smallScale.py > $@ 2> log/$@-err.txt) 2> log/$@-time.txt

main_table_allCNVs.txt: results.flag $D/mainTable-allCNVs.py
	(time $D/mainTable-allCNVs.py > $@ 2> log/$@-err.txt) 2> log/$@-time.txt

mainTable-simplex-multiplex-combined.txt: results.flag $D/mainTable-simplex-multiplex-combined.py
	(time $D/mainTable-simplex-multiplex-combined.py > $@ 2> log/$@-err.txt) 2> log/$@-time.txt

main_table_oneGeneCNVs.txt: results.flag $D/mainTable-oneGeneCNVs.py
	(time $D/mainTable-oneGeneCNVs.py > $@ 2> log/$@-err.txt) 2> log/$@-time.txt

main_table_nonCoding_%.txt: results.flag $D/mainTable.py
	(time $D/mainTable.py $* > $@ 2> log/$@-err.txt) 2> log/$@-time.txt

children_scatter.png: children_table.txt $D/figChildrenScatter.py
	(time $D/figChildrenScatter.py > log/$@-out.txt 2> log/$@-err.txt) 2> log/$@-time.txt

powerCNV.png: $D/figCNVPower.py
	(time $<  > log/$@-out.txt 2> log/$@-err.txt) 2> log/$@-time.txt

powerSnv.png: $D/figSnvPower.py ${DATA_FILES}
	(time $< $I > log/$@-out.txt 2> log/$@-err.txt) 2> log/$@-time.txt

CNVsByFilter.png: $D/figCNVsByFilter.py ${DATA_FILES}
	(time $<  > log/$@-out.txt 2> log/$@-err.txt) 2> log/$@-time.txt

figParentalAges.png: $D/figParentalAges.py ${DATA_FILES}
	(time $<  > log/$@-out.txt 2> log/$@-err.txt) 2> log/$@-time.txt

figRatesVsAge.png: $D/figRatesVsAge.py ${DATA_FILES}
	(time $<  > log/$@-out.txt 2> log/$@-err.txt) 2> log/$@-time.txt

ISB-power.png: $D/power-intronic-synonymous.py 
	(time $<  > log/$@-out.txt 2> log/$@-err.txt) 2> log/$@-time.txt

IRbyCNVGeneN.png: $D/figIRbyCNVgenes.py ${DATA_FILES}
	(time $<  > log/$@-out.txt 2> log/$@-err.txt) 2> log/$@-time.txt

ss-model-figures.flag: $D/splice_site_model.py ss-model-data.flag
	(time $<    > log/$@-out.txt 2> log/$@-err.txt && touch $@) 2> log/$@-time.txt

table_sanders_missing.txt: $D/table_sanders_missing.py ${DATA_FILES}
	(time $<    > log/$@-out.txt 2> log/$@-err.txt && touch $@) 2> log/$@-time.txt

propPNGs.flag: $D/draw_float_property_hist.py SSC_small_ann_funProp_denovo_table.txt
	(time (mkdir -p propPNGs && DI_DATA_EVENT_FILES=SSC_small_ann_funProp_denovo_table.txt $< propPNGs  > log/$@-out.txt 2> log/$@-err.txt && touch $@)) 2> log/$@-time.txt

######
######
## should probably remove
######
######
CDintron-list.txt: $D/list_introns.py
	(time $<    > $@ 2> log/$@-err.txt) 2> log/$@-time.txt

ss-model-data.flag: $D/prepare_splice_site_data.py  CDintron-list.txt
	(time $<    > log/$@-out.txt 2> log/$@-err.txt && touch $@) 2> log/$@-time.txt


