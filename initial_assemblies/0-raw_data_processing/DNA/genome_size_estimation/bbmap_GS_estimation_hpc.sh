#!/bin/bash
#
#SBATCH --job-name=R-BBMAP
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mail-user m.eitel@lmu.de
#SBATCH --mail-type=END,FAIL
#SBATCH --qos=low
#SBATCH --mem=50G


###############
### PRESETS ###
###############


R1=/home/meitel/data/avas/data/DNA/PE/Avastus_S0_L001_R1.fastq.gz
R2=/home/meitel/data/avas/data/DNA/PE/Avastus_S0_L001_R2.fastq.gz

NUMTH=12


#################
### BBmap run ###
#################



kmercountexact.sh \
in=$R1 \
in=$R2 \
k=31 \
out=Aphrocallistes_kmer_bbmap.kmer.counts \
threads=$NUMTH \
ploidy=4 \
-Xmx50G


