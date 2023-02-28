#!/bin/bash

################################################

#### Get average coverage of contigs after mapping the raw reads to it. To be use later for 
#### ploting Coverage vs GC% for extraction of bact contigs 

#### get_ave_cov.sh

################################################


# to print name, number of CG, number of total bases and GC% use this
# END { print h"\t"cg"\t"t"\t"(cg/t);}' filename

for x in depth.Aphrocallistes_reference_assembly_v1_ONT-RNA_gmap_renamed.bam ;do

    #name for final file
    a=$( echo "$x" | cut -f 2- -d "." ) 
    
    #extract names of contigs
    cat $x | awk '{print $1}' | sort | uniq >names.$a

    #Print columns names
    echo -e "Name_of_contig\tAverage_coverage\tSum_of_coverages\tNumber_of_bases" >>$a.ave.cov.gc.txt

    for y in `cat names.$a`;do

        #name of contig
        b=$( echo "$y" )
        
        #sum of coverage of the bases
        c=$( cat $x | grep $y | awk '{sum+=$3} END {print sum}' ) 

        #number of bases
        d=$( cat $x | grep $y | awk '{sum+=1} END {print sum}' ) 

        #average coverage
        e=$( echo "$c/$d" | bc )

        #print into a file
        #name of contig, average coverage, gc content, sum of coverage, number of bases
        echo -e "$b\t$e\t$c\t$d" >>$a.ave.cov.gc.txt

    done

done

#END
