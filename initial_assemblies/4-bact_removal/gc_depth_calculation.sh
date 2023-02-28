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

# to print name, number of CG, number of total bases and GC% use this
# END { print h"\t"cg"\t"t"\t"(cg/t);}' filename

for x in `ls -1 depth.*`;do

    #name for final file
    a=$( echo "$x" | cut -f 2- -d "." ) 
    
    #extract names of contigs
    cat $x | awk '{print $1}' | sort | uniq >names.$a

    #Print columns names
    echo -e "Name_of_contig\tAverage_coverage\tGC_content\tSum_of_coverages\tNumber_of_bases" >>$a.ave.cov.gc.txt

    for y in `cat names.$a`;do

        #name of contig
        b=$( echo "$y" )
        
        #sum of coverage of the bases
        c=$( cat $x | grep $y | awk '{sum+=$3} END {print sum}' ) 

        #number of bases
        d=$( cat $x | grep $y | awk '{sum+=1} END {print sum}' ) 

        #average coverage
        e=$( echo "$c/$d" | bc )

        #extract GC content of the contig from file
        f=$( cat GC.bact.txt | grep $y | awk '{print $2}' )

        #print into a file
        #name of contig, average coverage, gc content, sum of coverage, number of bases
        echo -e "$b\t$e\t$f\t$c\t$d" >>$a.ave.cov.gc.txt

    done

done

#END
