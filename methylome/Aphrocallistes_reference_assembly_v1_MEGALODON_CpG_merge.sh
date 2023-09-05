#!/bin/bash

NUMTH=40

DATASET=Aphrocallistes_reference_assembly_v1_MEGALODON_CpG
REF=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted.fasta

megalodon_extras merge modified_bases \
--max-processes $NUMTH \
--output-megalodon-results-dir Aphrocallistes_reference_assembly_v1_MEGALODON_CpG \
20171103_1624_MinION_megalodon_min_modbases_5mC_CpG_v001 \
20171109_1116_MinION_megalodon_min_modbases_5mC_CpG_v001 \
20180308_1413_PromethION_megalodon_min_modbases_5mC_CpG_v001 \
20200221_1309_MinION_megalodon_min_modbases_5mC_CpG_v001



## thr0.80
THRES1=0.80
megalodon_extras aggregate run --processes $NUMTH --megalodon-directory Aphrocallistes_reference_assembly_v1_MEGALODON_CpG --outputs mods --mod-output-formats bedmethyl --output-suffix thresh_${THRES1} --mod-binary-threshold $THRES1

## thr0.75
THRES2=0.75
megalodon_extras aggregate run --processes $NUMTH --megalodon-directory Aphrocallistes_reference_assembly_v1_MEGALODON_CpG --outputs mods --mod-output-formats bedmethyl --output-suffix thresh_${THRES2} --mod-binary-threshold $THRES2

## thr0.70
THRES3=0.70
megalodon_extras aggregate run --processes $NUMTH --megalodon-directory Aphrocallistes_reference_assembly_v1_MEGALODON_CpG --outputs mods --mod-output-formats bedmethyl --output-suffix thresh_${THRES3} --mod-binary-threshold $THRES3

## thr0.65
THRES4=0.65
megalodon_extras aggregate run --processes $NUMTH --megalodon-directory Aphrocallistes_reference_assembly_v1_MEGALODON_CpG --outputs mods --mod-output-formats bedmethyl --output-suffix thresh_${THRES4} --mod-binary-threshold $THRES4

## thr0.60
THRES5=0.60
megalodon_extras aggregate run --processes $NUMTH --megalodon-directory Aphrocallistes_reference_assembly_v1_MEGALODON_CpG --outputs mods --mod-output-formats bedmethyl --output-suffix thresh_${THRES5} --mod-binary-threshold $THRES5




# sorting of bed files
#grep ">" $REF | sed 's/>//' | cut -f 1 -d ' ' > ${REF%.fasta}.scaffold.ids

bedtools sort -faidx ${REF%.fasta}.scaffold.ids -i Aphrocallistes_reference_assembly_v1_MEGALODON_CpG/modified_bases.thresh_0.80.5mC.bed \
> ${DATASET}_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed

bedtools sort -faidx ${REF%.fasta}.scaffold.ids -i Aphrocallistes_reference_assembly_v1_MEGALODON_CpG/modified_bases.thresh_0.75.5mC.bed \
> ${DATASET}_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.75.5mC.bed

bedtools sort -faidx ${REF%.fasta}.scaffold.ids -i Aphrocallistes_reference_assembly_v1_MEGALODON_CpG/modified_bases.thresh_0.70.5mC.bed \
> ${DATASET}_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.70.5mC.bed

bedtools sort -faidx ${REF%.fasta}.scaffold.ids -i Aphrocallistes_reference_assembly_v1_MEGALODON_CpG/modified_bases.thresh_0.65.5mC.bed \
> ${DATASET}_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.65.5mC.bed

bedtools sort -faidx ${REF%.fasta}.scaffold.ids -i Aphrocallistes_reference_assembly_v1_MEGALODON_CpG/modified_bases.thresh_0.60.5mC.bed \
> ${DATASET}_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.60.5mC.bed

# cut percentage
cut -f 11 Aphrocallistes_reference_assembly_v1_MEGALODON_CpG/*0.80.5mC.bed | sort -n | uniq -c > Aphrocallistes_reference_assembly_v1_MEGALODON_CpG_v001_modified_bases.thresh_0.80.5mC.counts
cut -f 11 Aphrocallistes_reference_assembly_v1_MEGALODON_CpG/*0.75.5mC.bed | sort -n | uniq -c > Aphrocallistes_reference_assembly_v1_MEGALODON_CpG_v001_modified_bases.thresh_0.75.5mC.counts
cut -f 11 Aphrocallistes_reference_assembly_v1_MEGALODON_CpG/*0.70.5mC.bed | sort -n | uniq -c > Aphrocallistes_reference_assembly_v1_MEGALODON_CpG_v001_modified_bases.thresh_0.70.5mC.counts
cut -f 11 Aphrocallistes_reference_assembly_v1_MEGALODON_CpG/*0.65.5mC.bed | sort -n | uniq -c > Aphrocallistes_reference_assembly_v1_MEGALODON_CpG_v001_modified_bases.thresh_0.65.5mC.counts
cut -f 11 Aphrocallistes_reference_assembly_v1_MEGALODON_CpG/*0.60.5mC.bed | sort -n | uniq -c > Aphrocallistes_reference_assembly_v1_MEGALODON_CpG_v001_modified_bases.thresh_0.60.5mC.counts
