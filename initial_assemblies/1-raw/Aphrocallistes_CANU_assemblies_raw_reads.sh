#!/bin/bash
#
#SBATCH --job-name=CANU
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --clusters=mpp2

conda activate genomics

NUMTH=28

READS1=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_combined_5kb.fasta.gz
READS2=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_combined_10kb.fasta.gz
READS3=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_combined_15kb.fasta.gz

# CANU-1_5kb
canu -p Aphrocallistes_CANU_assembly-1_5kb -d Aphrocallistes_CANU_assembly-1_5kb \
corOutCoverage=40 \
ovlErrorRate=0.15 obtErrorRate=0.15 \
corMinCoverage=0 \
rawErrorRate=0.35 \
correctedErrorRate=0.15 \
minReadLength=5000 \
minOverlapLength=2000 \
corOverlapper=minimap \
genomeSize=80m \
maxMemory=350 \
maxThreads=$NUMTH \
useGrid=false \
-nanopore-raw $READS1 \
2> Aphrocallistes_CANU_assembly-1_5kb.log


# CANU-2_10kb
canu -p Aphrocallistes_CANU_assembly-2_10kb -d Aphrocallistes_CANU_assembly-2_10kb \
corOutCoverage=40 \
ovlErrorRate=0.15 obtErrorRate=0.15 \
corMinCoverage=0 \
rawErrorRate=0.35 \
correctedErrorRate=0.15 \
minReadLength=10000 \
minOverlapLength=2000 \
corOverlapper=minimap \
genomeSize=80m \
maxMemory=350 \
maxThreads=$NUMTH \
useGrid=false \
-nanopore-raw $READS2 \
2> Aphrocallistes_CANU_assembly-2_10kb.log


# CANU-3_15kb
canu -p Aphrocallistes_CANU_assembly-3_15kb -d Aphrocallistes_CANU_assembly-3_15kb \
corOutCoverage=40 \
ovlErrorRate=0.15 obtErrorRate=0.15 \
corMinCoverage=0 \
rawErrorRate=0.35 \
correctedErrorRate=0.15 \
minReadLength=15000 \
minOverlapLength=2000 \
corOverlapper=minimap \
stopOnLowCoverage=7 \
genomeSize=80m \
maxMemory=60 \
maxThreads=$NUMTH \
useGrid=false \
-nanopore-raw $READS3 \
2> Aphrocallistes_CANU_assembly-3_15kb.log


conda deactivate
