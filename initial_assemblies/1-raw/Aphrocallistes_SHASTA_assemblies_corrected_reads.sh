#!/bin/bash

source activate genomics

NUMTH=20

READS1=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_combined_5kb.fasta.gz
READS2=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_combined_10kb.fasta.gz
READS3=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_combined_15kb.fasta.gz

# SHASTA-1_5kb:
shastaGpu --command assemble \
--assemblyDirectory Aphrocallistes_SHASTA_assembly-1_5kb \
--memoryMode filesystem --memoryBacking 2M \
--Reads.minReadLength 5000 \
--Assembly.storeCoverageData \
--threads $NUMTH \
--useGpu \
--input $READS1

# SHASTA-2_10kb:
shastaGpu --command assemble \
--assemblyDirectory Aphrocallistes_SHASTA_assembly-2_10kb \
--memoryMode filesystem --memoryBacking 2M \
--Reads.minReadLength 10000 \
--Assembly.storeCoverageData \
--threads $NUMTH \
--useGpu \
--input $READS2

# SHASTA-2_15kb:
shastaGpu --command assemble \
--assemblyDirectory Aphrocallistes_SHASTA_assembly-2_15kb \
--memoryMode filesystem --memoryBacking 2M \
--Reads.minReadLength 15000 \
--Assembly.storeCoverageData \
--threads $NUMTH \
--useGpu \
--input $READS3

conda deactivate
