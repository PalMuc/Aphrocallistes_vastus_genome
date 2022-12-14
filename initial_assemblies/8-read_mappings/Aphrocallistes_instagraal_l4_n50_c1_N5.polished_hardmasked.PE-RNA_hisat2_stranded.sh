#!/bin/bash
#
#SBATCH --job-name=SRNA_HISAT2
#SBATCH --cluster=cm2_tiny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user m.eitel@lmu.de
#SBATCH --mail-type=END,FAIL


R1=/dss/dssfs02/lwp-dss-0001/uk213/uk213-dss-0000/di29vos/avas/data/RNA/PE/Avastus_100PE_RNAseq_stranded_f25_reheader_R1.fastq.gz
R2=/dss/dssfs02/lwp-dss-0001/uk213/uk213-dss-0000/di29vos/avas/data/RNA/PE/Avastus_100PE_RNAseq_stranded_f25_reheader_R2.fastq.gz

NUMTH=28



# short read mapping with HISAT2:
#for ASSEMBLY in *hardmasked.fasta ; do hisat2-build -p $NUMTH $ASSEMBLY ${ASSEMBLY%.fasta} ; done

for ASSEMBLY in *hardmasked.fasta ; \
do hisat2 -x ${ASSEMBLY%.fasta} -1 $R1 -2 $R2 -p $NUMTH --rna-strandness RF --max-intronlen 50000 --no-unal --no-repeat-index --no-mixed --no-discordant 2> ${ASSEMBLY%.fasta}.PE-RNA_hisat2.stranded.sorted.log |
samtools view -@ 12 -F 4 -Sbh | samtools sort -@ 12 -O BAM > ${ASSEMBLY%.fasta}.PE-RNA_hisat2_stranded.sorted.bam ; done
