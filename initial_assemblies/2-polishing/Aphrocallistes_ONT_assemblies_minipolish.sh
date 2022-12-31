#!/bin/bash


NUMTH=20

corREADS=/home/ubuntu/avas/data/DNA/ONT/Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_combined_5kb.fasta.gz
rawREADS=/home/ubuntu/avas/data/DNA/ONT/Aphrocallistes-DNA-guppy-3.2.4-R9.4.1-450bps_hac_porechopped_combined_5kb.fasta.gz

# corrected reads:
for GFA in *corrected_reads.gfa ; do minipolish-runner.py -t $NUMTH $corREADS $GFA > ${GFA%.gfa}.minipolished.gfa 2> ${GFA%.gfa}.minipolished.log ; done


# raw reads:
for GFA in *raw_reads.gfa ; do minipolish-runner.py -t $NUMTH $rawREADS $GFA > ${GFA%.gfa}.minipolished.gfa 2> ${GFA%.gfa}.minipolished.log ; done

