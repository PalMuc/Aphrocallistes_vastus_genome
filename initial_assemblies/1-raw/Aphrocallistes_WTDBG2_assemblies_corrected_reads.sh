#!/bin/bash
#
#SBATCH --job-name=WTDBG2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --clusters=mpp2

conda activate genomics

NUMTH=28

READS1=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_combined_5kb.fasta.gz
READS2=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_combined_10kb.fasta.gz
READS3=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_combined_15kb.fasta.gz

# WTDBG2-1_5kb:
## assemble long reads
wtdbg2 -x ont -g 80m -t $NUMTH \
-fo Aphrocallistes_WTDBG2_assembly-1_5kb \
-i $READS1

## derive consensus
wtpoa-cns -t NUMTH -i Aphrocallistes_WTDBG2_assembly-1_5kb.ctg.lay.gz -fo Aphrocallistes_WTDBG2_assembly-1_5kb.fasta


# WTDBG2-2_10kb:
## assemble long reads
wtdbg2 -x ont -g 80m -t $NUMTH \
-fo Aphrocallistes_WTDBG2_assembly-2_10kb \
-i $READS2

## derive consensus
wtpoa-cns -t NUMTH -i Aphrocallistes_WTDBG2_assembly-2_10kb.ctg.lay.gz -fo Aphrocallistes_WTDBG2_assembly-2_10kb.fasta


# WTDBG2-3_15kb:
## assemble long reads
wtdbg2 -x ont -g 80m -t $NUMTH \
-fo Aphrocallistes_WTDBG2_assembly-3_15kb \
-i $READS3

## derive consensus
wtpoa-cns -t NUMTH -i Aphrocallistes_WTDBG2_assembly-3_15kb.ctg.lay.gz -fo Aphrocallistes_WTDBG2_assembly-3_15kb.fasta



conda deactivate
