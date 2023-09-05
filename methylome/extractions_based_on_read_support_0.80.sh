#!/bin/bash


# 0.80 threshold results:
# extract positions that are supported by 100% of reads:
sort -k11,11 -n -r 20171103_1624_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 36133 | bedtools sort -i - > 20171103_1624_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.100.bedmethyl
sort -k11,11 -n -r 20171109_1116_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 31386 | bedtools sort -i - > 20171109_1116_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.100.bedmethyl
sort -k11,11 -n -r 20180308_1413_PromethION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 1887 | bedtools sort -i - > 20180308_1413_PromethION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.100.bedmethyl
sort -k11,11 -n -r 20200221_1309_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 19168 | bedtools sort -i - > 20200221_1309_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.100.bedmethyl
sort -k11,11 -n -r Aphrocallistes_MERGED_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 783 | bedtools sort -i - > Aphrocallistes_MERGED_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.100.bedmethyl

for file in *100.bedmethyl ; do cut -f 1-9 $file > filtered_beds/${file%.bedmethyl}.bed ; done


# extract positions that are supported by 50% up of reads:
sort -k11,11 -n -r 20171103_1624_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 106104 | bedtools sort -i - > 20171103_1624_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.50-up.bedmethyl
sort -k11,11 -n -r 20171109_1116_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 225219 | bedtools sort -i - > 20171109_1116_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.50-up.bedmethyl
sort -k11,11 -n -r 20180308_1413_PromethION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 23308 | bedtools sort -i - > 20180308_1413_PromethION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.50-up.bedmethyl
sort -k11,11 -n -r 20200221_1309_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 82668 | bedtools sort -i - > 20200221_1309_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.50-up.bedmethyl
sort -k11,11 -n -r Aphrocallistes_MERGED_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 16161 | bedtools sort -i - > Aphrocallistes_MERGED_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.50-up.bedmethyl

for file in *50-up.bedmethyl; do cut -f 1-9 $file > filtered_beds/${file%.bedmethyl}.bed ; done


# extract positions that are supported by 25% up of reads:
sort -k11,11 -n -r 20171103_1624_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 224912 | bedtools sort -i - > 20171103_1624_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.25-up.bedmethyl
sort -k11,11 -n -r 20171109_1116_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 575887 | bedtools sort -i - > 20171109_1116_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.25-up.bedmethyl
sort -k11,11 -n -r 20180308_1413_PromethION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 142715 | bedtools sort -i - > 20180308_1413_PromethION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.25-up.bedmethyl
sort -k11,11 -n -r 20200221_1309_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 230624 | bedtools sort -i - > 20200221_1309_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.25-up.bedmethyl
sort -k11,11 -n -r Aphrocallistes_MERGED_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 142666 | bedtools sort -i - > Aphrocallistes_MERGED_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.25-up.bedmethyl

for file in *25-up.bedmethyl; do cut -f 1-9 $file > filtered_beds/${file%.bedmethyl}.bed ; done


# extract positions that are supported by 20% up of reads:
sort -k11,11 -n -r 20171103_1624_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 260311 | bedtools sort -i - > 20171103_1624_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.20-up.bedmethyl
sort -k11,11 -n -r 20171109_1116_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 713182 | bedtools sort -i - > 20171109_1116_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.20-up.bedmethyl
sort -k11,11 -n -r 20180308_1413_PromethION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 214267 | bedtools sort -i - > 20180308_1413_PromethION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.20-up.bedmethyl
sort -k11,11 -n -r 20200221_1309_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 300440 | bedtools sort -i - > 20200221_1309_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.20-up.bedmethyl
sort -k11,11 -n -r Aphrocallistes_MERGED_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 225554 | bedtools sort -i - > Aphrocallistes_MERGED_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.20-up.bedmethyl

for file in *20-up.bedmethyl; do cut -f 1-9 $file > filtered_beds/${file%.bedmethyl}.bed ; done


# extract positions that are supported by 10% up of reads:
sort -k11,11 -n -r 20171103_1624_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 303961 | bedtools sort -i - > 20171103_1624_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.10-up.bedmethyl
sort -k11,11 -n -r 20171109_1116_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 1202901 | bedtools sort -i - > 20171109_1116_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.10-up.bedmethyl
sort -k11,11 -n -r 20180308_1413_PromethION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 559460 | bedtools sort -i - > 20180308_1413_PromethION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.10-up.bedmethyl
sort -k11,11 -n -r 20200221_1309_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 516867 | bedtools sort -i - > 20200221_1309_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.10-up.bedmethyl
sort -k11,11 -n -r Aphrocallistes_MERGED_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 642468 | bedtools sort -i - > Aphrocallistes_MERGED_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.10-up.bedmethyl

for file in *10-up.bedmethyl; do cut -f 1-9 $file > filtered_beds/${file%.bedmethyl}.bed ; done


# extract positions that are supported by 1% up of reads:
sort -k11,11 -n -r 20171103_1624_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 308141 | bedtools sort -i - > 20171103_1624_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.1-up.bedmethyl
sort -k11,11 -n -r 20171109_1116_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 1552037 | bedtools sort -i - > 20171109_1116_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.1-up.bedmethyl
sort -k11,11 -n -r 20180308_1413_PromethION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 1770952 | bedtools sort -i - > 20180308_1413_PromethION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.1-up.bedmethyl
sort -k11,11 -n -r 20200221_1309_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 558804 | bedtools sort -i - > 20200221_1309_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.1-up.bedmethyl
sort -k11,11 -n -r Aphrocallistes_MERGED_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.bed | head -n 2733082 | bedtools sort -i - > Aphrocallistes_MERGED_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.1-up.bedmethyl


for file in *1-up.bedmethyl; do cut -f 1-9 $file > filtered_beds/${file%.bedmethyl}.bed ; done
