#!/bin/bash

###############
### PRESETS ###
###############




WORKDIR=/home/ubuntu/avas/final_analyses/D_GENE_MODELLING/scaffolded_HiC
REF=${WORKDIR}/Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.fasta
BAM=${WORKDIR}/Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.sorted.bam
PROTEOME=${WORKDIR}/Aphrocallistes_RNAseq_Ind_1-5_combined_f25_reheader.combined.okay.fa.transdecoder.pep



NUMTH=20

export GENEMARK_PATH=/home/ubuntu/tools/braker2_dependencies/gmes_linux_64

###################################################################################
## long cRNA read hints & protein hints from the TransPi reference transcriptome ##
###################################################################################

## MINIMAP2 mapped reads:

braker.pl \
--workingdir=${WORKDIR}/01_BRAKER2/ONT-RNA/minimap2 \
--cores $NUMTH \
--AUGUSTUS_ab_initio \
--species=Aphrocallistes_instagraal_ONT-RNA_minimap2 \
--softmasking \
--gff3 \
--UTR=on \
--alternatives-from-evidence=true \
--augustus_args="--protein=on" \
--genome=$REF \
--bam=$BAM \
--prot_seq=$PROTEOME \
--prg=gth --gth2traingenes \
--makehub --email=m.eitel@lmu.de


