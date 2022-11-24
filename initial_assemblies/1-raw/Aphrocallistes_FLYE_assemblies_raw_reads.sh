#!/bin/bash
#
#SBATCH --job-name=FLYE
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --clusters=mpp2

conda activate genomics

NUMTH=28

READS1=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_combined_5kb.fasta.gz
READS2=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_combined_10kb.fasta.gz
READS3=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_combined_15kb.fasta.gz


# FLYE-1_5kb
flye --threads $NUMTH --genome-size 80m -i 2 \
--out-dir ~/Aphrocallistes_FLYE_assembly-1_5kb \
--nano-raw $READS1

# FLYE-2_10kb
flye --threads $NUMTH --genome-size 80m -i 2 \
--out-dir ~/Aphrocallistes_FLYE_assembly-2_10kb \
--nano-raw $READS2

# FLYE-3_15kb
flye --threads $NUMTH --genome-size 80m -i 2 \
--out-dir ~/Aphrocallistes_FLYE_assembly-3_15kb \
--nano-raw $READS3

conda deactivate
