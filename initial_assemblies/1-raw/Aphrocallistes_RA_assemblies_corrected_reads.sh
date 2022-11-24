#!/bin/bash
#
#SBATCH --job-name=RA
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --clusters=mpp3

conda activate genomics

NUMTH=64

READS1=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_combined_5kb.fasta.gz
READS2=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_combined_10kb.fasta.gz
READS3=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_combined_15kb.fasta.gz

# RA-1_5kb:
ra -t $NUMTH -x ont $READS1 \
> Aphrocallistes_RA_assembly-1_5kb.fasta \
2> Aphrocallistes_RA_assembly-1_5kb.log

# RA-2_10kb:
ra -t $NUMTH -x ont $READS2 \
> Aphrocallistes_RA_assembly-2_10kb.fasta \
2> Aphrocallistes_RA_assembly-2_10kb.log

# RA-3_15kb:
ra -t $NUMTH -x ont $READS3 \
> Aphrocallistes_RA_assembly-3_15kb.fasta \
2> Aphrocallistes_RA_assembly-3_15kb.log


conda deactivate

