#!/bin/bash
#
#SBATCH --job-name=SRNA_STAR
#SBATCH --cluster=cm2_tiny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user m.eitel@lmu.de
#SBATCH --mail-type=END,FAIL

NUMTH=28

READDIR=/dss/dssfs02/lwp-dss-0001/uk213/uk213-dss-0000/di29vos/avas/data/RNA/PE
REF=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.fasta


# stranded libraries:
IND8R1=$READDIR/IND-8_Av_100PE_RNAseq_stranded_f25_reheader_R1.fastq
IND8R2=$READDIR/IND-8_Av_100PE_RNAseq_stranded_f25_reheader_R2.fastq




# index reference genome:
#STAR --runMode genomeGenerate --runThreadN $NUMTH \
#--genomeFastaFiles $REF \
#--genomeDir ${REF%.fasta}

# separate libary read mapping:
STAR --runMode alignReads --runThreadN $NUMTH \
--genomeDir ${REF%.fasta} \
--readFilesIn $IND8R1 $IND8R2 \
--outSAMattrRGline ID:ind8 \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.PE-RNA_star_test_



mv Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.PE-RNA_star_test_Aligned.sortedByCoord.out.bam Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.PE-RNA_star_stranded_sorted.bam