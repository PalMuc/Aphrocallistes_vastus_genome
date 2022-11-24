#!/bin/bash
#
#SBATCH --job-name=MINIASM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --clusters=mpp2

conda activate genomics

NUMTH=28

READS1=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_combined_5kb.fasta.gz
READS2=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_combined_10kb.fasta.gz
READS3=Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_combined_15kb.fasta.gz


# MINIASM-1_5kb:
## all-against-all overlap:
minimap2 -t $NUMTH -x ava-ont \
$READS1 \
$READS1 \
2> Aphrocallistes_MINIASM_assembly-1_5kb_minimap2.log | gzip -1 > Aphrocallistes_MINIASM_assembly-1_5kb_minimap2.paf.gz

## assembly:
miniasm -f $READS1 \
Aphrocallistes_MINIASM_assembly-1_5kb_minimap2.paf.gz \
> Aphrocallistes_MINIASM_assembly-1_5kb.gfa \
2> Aphrocallistes_MINIASM_assembly-1_5kb.log

## convert GFA to fasta:
awk '/^S/{print ">"$2"\n"$3}' Aphrocallistes_MINIASM_assembly-1_5kb.gfa | fold > Aphrocallistes_MINIASM_assembly-1_5kb.fasta




# MINIASM-2_10kb:
## all-against-all overlap:
minimap2 -t $NUMTH -x ava-ont \
$READS2 \
$READS2 \
2> Aphrocallistes_MINIASM_assembly-2_10kb_minimap2.log | gzip -1 > Aphrocallistes_MINIASM_assembly-2_10kb_minimap2.paf.gz

## assembly:
miniasm -f $READS2 \
Aphrocallistes_MINIASM_assembly-2_10kb_minimap2.paf.gz \
> Aphrocallistes_MINIASM_assembly-2_10kb.gfa \
2> Aphrocallistes_MINIASM_assembly-2_10kb.log

## convert GFA to fasta:
awk '/^S/{print ">"$2"\n"$3}' Aphrocallistes_MINIASM_assembly-2_10kb.gfa | fold > Aphrocallistes_MINIASM_assembly-2_10kb.fasta




# MINIASM-3_15kb:
## all-against-all overlap:
minimap2 -t $NUMTH -x ava-ont \
$READS3 \
$READS3 \
2> Aphrocallistes_MINIASM_assembly-3_15kb_minimap2.log | gzip -1 > Aphrocallistes_MINIASM_assembly-3_15kb_minimap2.paf.gz

## assembly:
miniasm -f $READS3 \
Aphrocallistes_MINIASM_assembly-3_15kb_minimap2.paf.gz \
> Aphrocallistes_MINIASM_assembly-3_15kb.gfa \
2> Aphrocallistes_MINIASM_assembly-3_15kb.log

## convert GFA to fasta:
awk '/^S/{print ">"$2"\n"$3}' Aphrocallistes_MINIASM_assembly-3_15kb.gfa | fold > Aphrocallistes_MINIASM_assembly-3_15kb.fasta





conda deactivate
