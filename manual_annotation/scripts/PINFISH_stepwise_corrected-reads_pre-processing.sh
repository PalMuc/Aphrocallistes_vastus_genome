#!/bin/bash
#
#SBATCH --job-name=PIN_P
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --mem=100G
#SBATCH --mail-user=m.eitel@lmu.de
#SBATCH --mail-type=END,FAIL



###############
### PRESETS ###
###############



NUMTH=40

REF=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.fasta
CREADSA=/home/meitel/data/avas/data/RNA/ONT/lordec/Aphrocallistes_ONT_cDNA_run1-3_lordec_100bp.fasta
CREADSQ=Aphrocallistes_ONT_cDNA_run1-3_lordec_100bp.fastq
TEMP=/home/meitel/data/avas/pinfish/tmp

####################################
### data preparation for PINFISH ###
####################################

## Using BBmap to create false quality scores (before mapping).
reformat.sh in=$CREADSA out=$CREADSQ qfake=35

# map LORDEC corrected reads to the hardmasked genome with MINIMAP2:
minimap2 -ax splice -ub --secondary=no -t $NUMTH $REF $CREADSQ > ${REF%.fasta}.ONT-RNA_minimap2.sam


## Filtering of unmapped and secondary alignments.
## This script only filters reads found twice, it there are reads mapped more times they will affect the PINFISH run (Still needs to be improved to cover all secondary alignments).
## The -q 1 option on samtools view is not enough to filter, given that there are some duplicates that have quality score 1 and they will be kept.

# remove unmapped reads
samtools view -t $NUMTH -F 4 -h -S \
${REF%.fasta}.ONT-RNA_minimap2.sam > ${REF%.fasta}.ONT-RNA_minimap2_mapped.sam

# get the ids of reads with secondary alignments (duplicates)
samtools view -t $NUMTH -S \
${REF%.fasta}.ONT-RNA_minimap2_mapped.sam \
| cut -f 1 | sort | uniq -c | grep -v $'\t1 ' | awk '{print "\t2 "$0}' - > duplicates.ids

# remove duplicates
grep -Ff duplicates.ids -v \
${REF%.fasta}.ONT-RNA_minimap2_mapped.sam  \
> ${REF%.fasta}.ONT-RNA_minimap2_mapped_singletons_filtered.sam

# generate sorted bam file
samtools view -t $NUMTH -b \
${REF%.fasta}.ONT-RNA_minimap2_mapped_singletons_filtered.sam | samtools sort -t $NUMTH \
> ${REF%.fasta}.ONT-RNA_minimap2_mapped_singletons_filtered_sorted.bam

# check if all duplicates have been removed successfully
samtools view -t $NUMTH  \
${REF%.fasta}.ONT-RNA_minimap2_mapped_singletons_filtered_sorted.bam | cut -f 1 | sort | uniq -c | grep -c $'\t2 ' \
> duplicates.check

