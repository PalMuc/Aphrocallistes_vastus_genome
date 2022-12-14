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


# unstranded libraries:
IND1BR1=$READDIR/IND-1_Av2-3-B-2_S7_combined_f25_R1.fastq.gz
IND1BR2=$READDIR/IND-1_Av2-3-B-2_S7_combined_f25_R2.fastq.gz
IND1OR1=$READDIR/IND-1_AV-20-10_S1_combined_f25_R1.fastq.gz
IND1OR2=$READDIR/IND-1_AV-20-10_S1_combined_f25_R1.fastq.gz
IND2OR1=$READDIR/IND-2_AV-4-0-6_S12_combined_f25_R1.fastq.gz
IND2OR2=$READDIR/IND-2_AV-4-0-6_S12_combined_f25_R2.fastq.gz
IND3OR1=$READDIR/IND-3_AV3-0-6_S9_combined_f25_R1.fastq.gz
IND3OR2=$READDIR/IND-3_AV3-0-6_S9_combined_f25_R2.fastq.gz
IND4BR1=$READDIR/IND-4_B1-4-2_S10_combined_f25_R1.fastq.gz
IND4BR2=$READDIR/IND-4_B1-4-2_S10_combined_f25_R2.fastq.gz
IND5BR1=$READDIR/IND-5_B2-8_S5_combined_f25_R1.fastq.gz
IND5BR2=$READDIR/IND-5_B2-8_S5_combined_f25_R2.fastq.gz
IND6BR1=$READDIR/IND-6_AV4-1_R1.f25.fastq.gz
IND6BR2=$READDIR/IND-6_AV4-1_R2.f25.fastq.gz
IND7BR1=$READDIR/IND-7_Body_SRR1068281_t.f.25_R1.fastq.gz
IND7BR2=$READDIR/IND-7_Body_SRR1068281_t.f.25_R2.fastq.gz

# stranded libraries:
IND8BR1=$READDIR/IND-8_Av_100PE_RNAseq_stranded_f25_reheader_R1.fastq.gz
IND8BR2=$READDIR/IND-8_Av_100PE_RNAseq_stranded_f25_reheader_R2.fastq.gz




# index reference genome:
#STAR --runMode genomeGenerate --runThreadN $NUMTH \
#--genomeFastaFiles $REF \
#--genomeDir ${REF%.fasta}

# separate libary read mapping:
STAR --runMode alignReads --runThreadN $NUMTH \
--genomeDir ${REF%.fasta} \
--readFilesIn $INDB1R1,$INDO1R1,$INDO2R1,$IND3OR1,$IND4BR1,$IND5BR1,$IND6BR1,$IND7BR1,$IND8BR1 $INDB1R2,$INDO1R2,$INDO2R2,$IND3OR2,$IND4BR2,$IND5BR2,$IND6BR2,$IND7BR2,$IND8BR2 \
--readFilesCommand gunzip -c \
--outSAMattrRGline ID:ind1_Body , ID:ind1_Osculum , ID:ind2_Osculum , ID:ind3_Osculum , ID:ind4_Body , ID:ind5_Body , ID:ind6_Body , ID:ind7_Body , ID:ind8_Body \
--outSAMtype SAM \
--outFileNamePrefix Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.PE-RNA_star_

#convert sam to bam and sort
samtools view -F 4 -@ $NUMTH -b -S ${REF%.fasta}.PE-RNA_star_Aligned.out.sam | \
samtools sort -@ $NUMTH > ${REF%.fasta}..PE-RNA_star.sorted.bam ; done

