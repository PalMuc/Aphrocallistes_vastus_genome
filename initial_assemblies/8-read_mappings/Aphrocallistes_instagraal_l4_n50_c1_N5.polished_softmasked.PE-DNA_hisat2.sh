#!/bin/bash
#
#SBATCH --job-name=SDNA_HISAT2
#SBATCH --cluster=cm2_tiny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user m.eitel@lmu.de
#SBATCH --mail-type=END,FAIL


R1=/dss/dssfs02/lwp-dss-0001/uk213/uk213-dss-0000/di29vos/avas/data/DNA/PE/karect_Avastus_S0_L001_R1.t.f.fastq.gz
R2=/dss/dssfs02/lwp-dss-0001/uk213/uk213-dss-0000/di29vos/avas/data/DNA/PE/karect_Avastus_S0_L001_R2.t.f.fastq.gz

NUMTH=28



# short read mapping with HISAT2:
for ASSEMBLY in *softmasked.fasta ; do hisat2-build -p $NUMTH $ASSEMBLY ${ASSEMBLY%.fasta} ; done

for ASSEMBLY in *softmasked.fasta ; \
do hisat2 -x ${ASSEMBLY%.fasta} -1 $R1 -2 $R2 -p $NUMTH --no-spliced-alignment --no-unal --no-repeat-index --no-mixed --no-discordant 2> ${ASSEMBLY%.fasta}.PE-DNA_hisat2.log |
samtools view -@ 12 -F 4 -Sbh | samtools sort -@ 12 -O BAM > ${ASSEMBLY%.fasta}.PE-DNA_hisat2.bam ; done
