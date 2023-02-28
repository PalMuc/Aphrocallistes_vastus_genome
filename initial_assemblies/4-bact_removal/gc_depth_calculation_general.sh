#!/bin/bash

################################################

#### Get average coverage of contigs after mapping the raw reads to it. To be use later for 
#### ploting Coverage vs GC% for extraction of bact contigs 

#### get_ave_cov.sh

################################################

#BEGIN

#PATH to contigs file
contigs=/Volumes/Michis_Data/WORK/Genomes/AVAS/13_ANALYSES/A_CONTAMINATION_REMOVAL/Aphrocallistes_reference_assembly_v1.fasta
#PATH to bam files folder
bamf=/Volumes/Michis_Data/WORK/Genomes/AVAS/09_MAPPINGS/reference/

#Get names from bam files (all before the first .)
for x in $bamf/*.bam;do

	#name for the depth files
    a=$( echo "$x" | sed "s|$bamf\/||" | cut -f 1-6 -d "_" ) 
    
    #create a separate file for the depth results
    samtools depth $x >depth.$a 

done

#Calculate GC content
#To be use later for the loop
awk 'BEGIN { FS=""; cg=0; t=0; } \
{ \
    if ($1 != ">") { \
        for (i = 1; i <= NF; i++) { \
            if ($i ~ /[ACTGactg]/) { t++; } \
            if ($i ~ /[CGcg]/) { cg++; }}} \
    else { \
        if (t > 0) { \
            print h"\t"(cg/t); cg = 0; t = 0; } \
        h = substr($0,2); } \
} END { print h"\t"(cg/t);}' $contigs >GC.bact.txt

