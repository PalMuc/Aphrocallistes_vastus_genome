#!/bin/bash

NUMTH=10


DATASET=20171103_1624_MinION
FAST5=/home/ubuntu/avas/data/DNA/ONT/raw/MinION/Aphrocallistes_${DATASET}/fast5
REF=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted.fasta

megalodon \
--guppy-server-path=/home/ubuntu/tools/ont-guppy/bin/guppy_basecall_server \
--guppy-params "--num_callers 8 --ipc_threads 12 --chunk_size 100" \
--guppy-config res_dna_r941_min_modbases_5mC_CpG_v001.cfg \
--devices 0 1 \
--processes $NUMTH \
--outputs "mods" \
--mod-output-formats "bedmethyl" \
--mod-binary-threshold 0.8 \
--guppy-timeout 10.0 \
--output-directory ${DATASET}_megalodon_min_modbases_5mC_CpG_v001 \
--reference $REF \
--overwrite \
$FAST5



mv ${DATASET}_megalodon_min_modbases_5mC_CpG_v001/modified_bases.5mC.bed ${DATASET}_megalodon_min_modbases_5mC_CpG_v001/modified_bases.thresh_0.80.5mC.bed


#re-aggregate calls:

## thr0.75
THRES1=0.75
megalodon_extras aggregate run --processes $NUMTH --megalodon-directory ${DATASET}_megalodon_min_modbases_5mC_CpG_v001 --outputs mods --mod-output-formats bedmethyl --output-suffix thresh_${THRES1} --mod-binary-threshold $THRES1

## thr0.70
THRES2=0.70
megalodon_extras aggregate run --processes $NUMTH --megalodon-directory ${DATASET}_megalodon_min_modbases_5mC_CpG_v001 --outputs mods --mod-output-formats bedmethyl --output-suffix thresh_${THRES2} --mod-binary-threshold $THRES2

## thr0.65
THRES3=0.65
megalodon_extras aggregate run --processes $NUMTH --megalodon-directory ${DATASET}_megalodon_min_modbases_5mC_CpG_v001 --outputs mods --mod-output-formats bedmethyl --output-suffix thresh_${THRES3} --mod-binary-threshold $THRES3

## thr0.60
THRES4=0.60
megalodon_extras aggregate run --processes $NUMTH --megalodon-directory ${DATASET}_megalodon_min_modbases_5mC_CpG_v001 --outputs mods --mod-output-formats bedmethyl --output-suffix thresh_${THRES4} --mod-binary-threshold $THRES4




# sorting of bed files
grep ">" $REF | sed 's/>//' | cut -f 1 -d ' ' > ${REF%.fasta}.scaffold.ids

bedtools sort -faidx ${REF%.fasta}.scaffold.ids -i ${DATASET}_megalodon_min_modbases_5mC_CpG_v001/modified_bases.thresh_0.80.5mC.bed \
> ${DATASET}_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed

bedtools sort -faidx ${REF%.fasta}.scaffold.ids -i ${DATASET}_megalodon_min_modbases_5mC_CpG_v001/modified_bases.thresh_0.75.5mC.bed \
> ${DATASET}_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.75.5mC.bed

bedtools sort -faidx ${REF%.fasta}.scaffold.ids -i ${DATASET}_megalodon_min_modbases_5mC_CpG_v001/modified_bases.thresh_0.70.5mC.bed \
> ${DATASET}_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.70.5mC.bed

bedtools sort -faidx ${REF%.fasta}.scaffold.ids -i ${DATASET}_megalodon_min_modbases_5mC_CpG_v001/modified_bases.thresh_0.65.5mC.bed \
> ${DATASET}_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.65.5mC.bed

bedtools sort -faidx ${REF%.fasta}.scaffold.ids -i ${DATASET}_megalodon_min_modbases_5mC_CpG_v001/modified_bases.thresh_0.60.5mC.bed \
> ${DATASET}_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.60.5mC.bed

# cut percentage
cut -f 11 ${DATASET}_megalodon_min_modbases_5mC_CpG_v001/*0.80.5mC.bed | sort -n | uniq -c > ${DATASET}_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.counts
cut -f 11 ${DATASET}_megalodon_min_modbases_5mC_CpG_v001/*0.75.5mC.bed | sort -n | uniq -c > ${DATASET}_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.75.5mC.counts
cut -f 11 ${DATASET}_megalodon_min_modbases_5mC_CpG_v001/*0.70.5mC.bed | sort -n | uniq -c > ${DATASET}_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.70.5mC.counts
cut -f 11 ${DATASET}_megalodon_min_modbases_5mC_CpG_v001/*0.65.5mC.bed | sort -n | uniq -c > ${DATASET}_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.65.5mC.counts
cut -f 11 ${DATASET}_megalodon_min_modbases_5mC_CpG_v001/*0.60.5mC.bed | sort -n | uniq -c > ${DATASET}_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.60.5mC.counts

