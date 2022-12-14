#!/bin/bash

NUMTH=40
READS=Aphrocallistes_cDNA_combined_pychoppper_100bp_to_20kb.fasta.gz
ASSEMBLY=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.fasta

# to generate sam file (mapping of the dna genome)
minimap2 -ax splice -ub --secondary=no -t $NUMTH $ASSEMBLY $READS > ${ASSEMBLY%.fasta}.ONT-RNA_minimap2_raw-reads.sam

# convert sam to bam and sort
samtools view -h -F4 -@ $NUMTH -bS ${ASSEMBLY%.fasta}.ONT-RNA_minimap2_raw-reads.sam | samtools sort -@ $NUMTH > ${ASSEMBLY%.fasta}.ONT-RNA_minimap2_raw-reads.sorted.bam



