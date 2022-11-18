#!/bin/bash

PWD=/Users/meitel/Desktop/makehomologs/holozoa_allprots_v10
ANNTAB=$PWD/ids/avas_v1_29_annotations_prot_fr_blastp_k10.annotation_table.tsv
ANNTABHEADER=$PWD/ids/avas_v1_29_annotations_prot_fr_blastp_k10.annotation_table_header.tsv
HOMOMETALL=$PWD/ids/homologs_metazoa_all_clusters.txt
HOMOMETSOME=$PWD/ids/homologs_metazoa_some_clusters.txt
HOMOPORALL=$PWD/ids/homologs_porifera_all_clusters.txt
HOMOPORSOME=$PWD/ids/homologs_porifera_some_clusters.txt
HOMOHEXALL=$PWD/ids/homologs_hexactinellida_all_clusters.txt
HOMOHEXSOME=$PWD/ids/homologs_hexactinellida_some_clusters.txt
PARAMETALL=$PWD/ids/paralogs_metazoa_all_clusters.txt
PARAMETSOME=$PWD/ids/paralogs_metazoa_some_clusters.txt
PARAPORALL=$PWD/ids/paralogs_porifera_all_clusters.txt
PARAPORSOME=$PWD/ids/paralogs_porifera_some_clusters.txt
PARAHEXALL=$PWD/ids/paralogs_hexactinellida_all_clusters.txt
PARAHEXSOME=$PWD/ids/paralogs_hexactinellida_some_clusters.txt

mkdir $PWD/ids/homologs
mkdir $PWD/ids/paralogs
mkdir $PWD/ids/annotations_per_cluster
mkdir $PWD/ids/annotations_per_cluster/homologs
mkdir $PWD/ids/annotations_per_cluster/paralogs
mkdir $PWD/ids/annotations_per_cluster_final
mkdir $PWD/ids/annotations_per_cluster_final/homologs
mkdir $PWD/ids/annotations_per_cluster_final/paralogs

# get AVAS seq identities per cluster:
cd $PWD/clusters_holozoa_allprots_v10/homologs
for file in *.fasta ; do ls $file | awk '{print("mv "$1" "$1)}' | sed 's/homologs_holozoa_allprots_v10/homolog_cluster/2' | /bin/sh ; done
for file in *.fasta ; do grep '>' $file | sed 's/>//' | grep 'AVAS' | sed 's/AVAS_//' > ../../ids/homologs/${file%.fasta}.ids ; find ../../ids/homologs -empty -type f -delete ; done

cd $PWD/clusters_holozoa_allprots_v10/paralogs
for file in *.fasta ; do ls $file | awk '{print("mv "$1" "$1)}' | sed 's/paralogs_holozoa_allprots_v10/paralog_cluster/2' | /bin/sh ; done
for file in *.fasta ; do grep '>' $file | sed 's/>//' | grep 'AVAS' | sed 's/AVAS_//' > ../../ids/paralogs/${file%.fasta}.ids ; find ../../ids/paralogs -empty -type f -delete ; done



# get AVAS cluster lits for homologs:
cd $PWD/ids/homologs
ls | sed 's/\.ids//' > ../cluster_list_homologs.txt
for file in *ids ; do tr '\n' ';' < $file | sed 's/;$//' > ${file%.ids}_revised.ids ; done 
cat *_revised.ids > ../ids_list_homologs.txt

paste cluster_list_homologs.txt ids_list_homologs.txt > holozoa_allprots_v10_homologs.tsv

grep -Ff $HOMOMETALL holozoa_allprots_v10_homologs.tsv > holozoa_allprots_v10_homologs_metazoa_all.tsv
grep -Ff $HOMOMETSOME holozoa_allprots_v10_homologs.tsv > holozoa_allprots_v10_homologs_metazoa_some.tsv
grep -Ff $HOMOPORALL holozoa_allprots_v10_homologs.tsv > holozoa_allprots_v10_homologs_porifera_all.tsv
grep -Ff $HOMOPORSOME holozoa_allprots_v10_homologs.tsv > holozoa_allprots_v10_homologs_porifera_some.tsv
grep -Ff $HOMOHEXALL holozoa_allprots_v10_homologs.tsv > holozoa_allprots_v10_homologs_hexactinellida_all.tsv
grep -Ff $HOMOHEXSOME holozoa_allprots_v10_homologs.tsv > holozoa_allprots_v10_homologs_hexactinellida_some.tsv



# get AVAS cluster lits for paralogs:
cd $PWD/ids/paralogs
ls | sed 's/\.ids//' > ../cluster_list_paralogs.txt
for file in *ids ; do tr '\n' ';' < $file | sed 's/;$//' > ${file%.ids}_revised.ids ; done
cat *_revised.ids > ../ids_list_paralogs.txt

paste cluster_list_paralogs.txt ids_list_paralogs.txt > holozoa_allprots_v10_paralogs.tsv

grep -Ff $PARAMETALL holozoa_allprots_v10_paralogs.tsv > holozoa_allprots_v10_paralogs_metazoa_all.tsv
grep -Ff $PARAMETSOME holozoa_allprots_v10_paralogs.tsv > holozoa_allprots_v10_paralogs_metazoa_some.tsv
grep -Ff $PARAPORALL holozoa_allprots_v10_paralogs.tsv > holozoa_allprots_v10_paralogs_porifera_all.tsv
grep -Ff $PARAPORSOME holozoa_allprots_v10_paralogs.tsv > holozoa_allprots_v10_paralogs_porifera_some.tsv
grep -Ff $PARAHEXALL holozoa_allprots_v10_paralogs.tsv > holozoa_allprots_v10_paralogs_hexactinellida_all.tsv
grep -Ff $PARAHEXSOME holozoa_allprots_v10_paralogs.tsv > holozoa_allprots_v10_paralogs_hexactinellida_some.tsv



# extract functional annotation data and add to AVAS clusters:
cd $PWD/ids/homologs
for file in *ids ; do grep -Ff $file $ANNTAB > ../annotations_per_cluster/homologs/${file%.ids}.annot.tsv ; cat $ANNTABHEADER ../annotations_per_cluster/homologs/${file%.ids}.annot.tsv > ../annotations_per_cluster_final/homologs/${file%.ids}.annot.final.tsv ; done

cd $PWD/ids/paralogs
for file in *ids ; do grep -Ff $file $ANNTAB > ../annotations_per_cluster/paralogs/${file%.ids}.annot.tsv ; cat $ANNTABHEADER ../annotations_per_cluster/paralogs/${file%.ids}.annot.tsv > ../annotations_per_cluster_final/paralogs/${file%.ids}.annot.final.tsv ; done




