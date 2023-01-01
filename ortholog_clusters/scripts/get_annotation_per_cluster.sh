#!/bin/bash

PWD=/Volumes/Michis_Data/WORK/Genomes/AVAS/13_FINAL_ANALYSES/orthology_assignment/makehomologs/holozoa_allprots_v10
ANNTAB=$PWD/ids/avas_v1_29_annotations_prot_fr_blastp_k10.annotation_table.tsv
ANNTABHEADER=$PWD/ids/avas_v1_29_annotations_prot_fr_blastp_k10.annotation_table_header.tsv


# get AVAS seq identities per cluster:
mkdir $PWD/ids/paralogs
mkdir $PWD/ids/homologs

cd $PWD/clusters_holozoa_allprots_v10/paralogs
for file in *fasta ; do grep '>' $file | sed 's/>//' | grep 'AVAS' | sed 's/AVAS_//' > $PWD/ids/paralogs/avas.${file%.fasta.aln}.ids ; find $PWD/ids/paralogs/ -empty -type f -delete ; done

cd $PWD/clusters_holozoa_allprots_v10/homologs
for file in *fasta ; do grep '>' $file | sed 's/>//' | grep 'AVAS' | sed 's/AVAS_//' > $PWD/ids/homologs/avas.${file%.fasta.aln}.ids ; find $PWD/ids/homologs/ -empty -type f -delete ; done

# extract functional annotation data and add to clusters:
cd $PWD/ids/homologs
for file in *ids ; do grep -Ff $file $ANNTAB > ${file%.ids}.annot.tsv ; cat $ANNTABHEADER ${file%.ids}.annot.tsv > ${file%.tsv}.final.tsv ; done

cd $PWD/ids/paralogs
for file in *ids ; do grep -Ff $file $ANNTAB > ${file%.ids}.annot.tsv ; cat $ANNTABHEADER ${file%.ids}.annot.tsv > ${file%.tsv}.final.tsv ; done

mkdir $PWD/ids/annotations_per_cluster
mkdir $PWD/ids/annotations_per_cluster/homologs
mkdir $PWD/ids/annotations_per_cluster/paralogs
mkdir $PWD/ids/annotations_per_cluster_final
mkdir $PWD/ids/annotations_per_cluster_final/homologs
mkdir $PWD/ids/annotations_per_cluster_final/paralogs

for file in $PWD/ids/homologs/*.annot.tsv ; do mv $file $PWD/ids/annotations_per_cluster/homologs ; done
for file in $PWD/ids/homologs/*.final.tsv ; do mv $file $PWD/ids/annotations_per_cluster_final/homologs ; done
for file in $PWD/ids/paralogs/*.annot.tsv ; do mv $file $PWD/ids/annotations_per_cluster/paralogs ; done
for file in $PWD/ids/paralogs/*.final.tsv ; do mv $file $PWD/ids/annotations_per_cluster_final/paralogs ; done



