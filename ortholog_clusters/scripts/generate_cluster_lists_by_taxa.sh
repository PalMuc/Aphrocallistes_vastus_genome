#!/bin/bash

PWD=/Volumes/Michis_Data/WORK/Genomes/AVAS/13_FINAL_ANALYSES/orthology_assignment/makehomologs/holozoa_allprots_v10/ids/paralogs

for file in $PWD/*.ids ; do ls $file | awk '{print("mv "$1" "$1)}' | sed 's/avas.paralogs_holozoa_allprots_v10/paralogs_cluster/2' | /bin/sh ; done
for file in $PWD/*.ids ; do ls $file | awk '{print("mv "$1" "$1)}' | sed 's/.fasta//2' | /bin/sh ; done

ls *ids > cluster_list.txt
for file in $PWD/*ids ; do paste -d '\n' $file | sed -e 's/,$//' >> ids_list.txt ; done




paste cluster_list.txt ids_list.txt > holozoa_allprots_v10_clusters_porifera.tsv
paste cluster_list.txt ids_list.txt > holozoa_allprots_v10_clusters_metazoa.tsv
paste cluster_list.txt ids_list.txt > holozoa_allprots_v10_clusters_hexactinellida.tsv