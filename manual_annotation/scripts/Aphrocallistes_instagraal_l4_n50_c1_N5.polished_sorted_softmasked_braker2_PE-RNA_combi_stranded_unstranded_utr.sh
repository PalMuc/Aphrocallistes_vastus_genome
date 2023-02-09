#!/bin/bash

###############
### PRESETS ###
###############



WORKDIR=/home/ubuntu/avas/final_analyses/D_GENE_MODELLING/scaffolded_HiC
REF=${WORKDIR}/Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.fasta
FBAM=${WORKDIR}/Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.PE-RNA_hisat2_stranded.fwd.bam
RBAM=${WORKDIR}/Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.PE-RNA_hisat2_stranded.rev.bam
UBAM=${WORKDIR}/Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.PE-RNA_hisat2.sorted.bam
PROTEOME=${WORKDIR}/Aphrocallistes_RNAseq_Ind_1-5_combined_f25_reheader.combined.okay.fa.transdecoder.pep

NUMTH=20

export GENEMARK_PATH=/home/ubuntu/tools/braker2_dependencies/gmes_linux_64

#######################################################################################
## short RNA-seq read hints & protein hints from the TransPi reference transcriptome ##
#######################################################################################



braker.pl \
--workingdir=${WORKDIR}/01_BRAKER2/PE-RNA_combi_stranded_unstranded_utr \
--cores $NUMTH \
--AUGUSTUS_ab_initio \
--species=Aphrocallistes_instagraal_PE-RNA_combi_stranded_unstranded_utr \
--softmasking \
--gff3 \
--UTR=on \
--augustus_args="--protein=on" \
--genome=$REF \
--bam=$FBAM,$RBAM,$UBAM \
--stranded=+,-,. \
--prot_seq=$PROTEOME \
--prg=gth --gth2traingenes \
--makehub --email=m.eitel@lmu.de
