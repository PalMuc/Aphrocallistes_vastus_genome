#!/bin/bash
#
#SBATCH --job-name=SDNA_BOWTIE2
#SBATCH --cluster=cm2_tiny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user m.eitel@lmu.de
#SBATCH --mail-type=END,FAIL


R1=/dss/dssfs02/lwp-dss-0001/uk213/uk213-dss-0000/di29vos/avas/data/DNA/PE/karect_Avastus_S0_L001_R1.t.f.fastq.gz
R2=/dss/dssfs02/lwp-dss-0001/uk213/uk213-dss-0000/di29vos/avas/data/DNA/PE/karect_Avastus_S0_L001_R2.t.f.fastq.gz

NUMTH=28


# build index Files:
for ASSEMBLY in *softmasked.fasta ; do bowtie2-build --threads $NUMTH $ASSEMBLY ${ASSEMBLY%.fasta} ; done

#  mapping corrected PE-reads:
for ASSEMBLY in *softmasked.fasta ; do \
bowtie2 -x ${ASSEMBLY%.fasta} -1 $R1 -2 $R2 --threads $NUMTH --sensitive-local --end-to-end --no-mixed --no-discordant --no-unal \
2> ${ASSEMBLY%.fasta}.PE-DNA_bowtie2.log > ${ASSEMBLY%.fasta}.PE-DNA_bowtie2.sam ; done

#convert sam to bam and sort
for ASSEMBLY in *softmasked.fasta ; do samtools view -F 4 -@ 12 -b -S ${ASSEMBLY%.fasta}.PE-DNA_bowtie2.sam | \
samtools sort -@ 12 > ${ASSEMBLY%.fasta}.PE-DNA_bowtie2.sorted.bam ; done

# rm *sam


