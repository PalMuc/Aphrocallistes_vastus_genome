#!/bin/bash

REF=Aphrocallistes_instagraal_l4_n50_c1_N5.polished.fasta

# sort by length, rename, remove starting/ending "N", remove scaffolds ,1kb:

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $REF | \
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |sort -k1,1nr | cut -f 2- | tr "\t" "\n" | awk '/^>/{print ">" ++i; next}{print}' | \
fastarenamer.py -p "Aphrocalllistes_vastus_HiC-scaffold" - | \
sed -E 's/^[N]+//' | sed -E 's/[N]+$//' | \
sed -E 's/Aphrocalllistes_vastus_HiC-scaffold_([0-9]{2}$)/Aphrocalllistes_vastus_HiC-scaffold_0\1/' | \
sed -E 's/Aphrocalllistes_vastus_HiC-scaffold_([0-9]{1}$)/Aphrocalllistes_vastus_HiC-scaffold_00\1/' | \
sizecutter.py -n -a 1000 -p - 2> ${REF%.fasta}_sorted.fasta.stats > ${REF%.fasta}_sorted.fasta
