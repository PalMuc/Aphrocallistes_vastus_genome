#!/bin/bash
#
#SBATCH --job-name=RAVEN
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --clusters=mpp2

conda activate genomics

NUMTH=28

READS1=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_combined_5kb.fasta.gz
READS2=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_combined_10kb.fasta.gz
READS3=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_combined_15kb.fasta.gz

# RAVEN-1_5kb:
raven -t $NUMTH -p 0 $READS1
> Aphrocallistes_RAVEN_assembly-1_5kb.fasta \
2> Aphrocallistes_RAVEN_assembly-1_5kb.log

# RAVEN-2_10kb:
raven -t $NUMTH -p 0 $READS2
> Aphrocallistes_RAVEN_assembly-2_10kb.fasta \
2> Aphrocallistes_RAVEN_assembly-2_10kb.log

# RAVEN-3_15kb:
raven -t $NUMTH -p 0 $READS3
> Aphrocallistes_RAVEN_assembly-3_15kb.fasta \
2> Aphrocallistes_RAVEN_assembly-3_15kb.log

conda deactivate