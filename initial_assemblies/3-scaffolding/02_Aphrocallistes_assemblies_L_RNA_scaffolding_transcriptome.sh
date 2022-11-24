#!/bin/bash


# https://github.com/CAFS-bioinformatics/L_RNA_scaffolder


READS=Aphrocallistes_RNAseq_Ind_1-5_combined_f25_reheader.combined.okay.fasta

NUMTH=10

source activate genomics

# generate .psl BLAT alignment file:
for POLISHED in *fasta ; do pblat -t=dna -q=dna -noHead -minIdentity=95 -maxIntron=50000 -threads=$NUMTH $POLISHED $READS ${POLISHED%.fasta}.BLAT_ONT.psl ; echo "done" ; done

conda deactivate

# running L_RNA_scaffolder
for POLISHED in *.fasta ; do mkdir ${POLISHED%.fasta}_scaffolded ; \
do /home/ubuntu/tools/L_RNA_scaffolder/L_RNA_scaffolder.sh -f 1 -o ${POLISHED%.fasta}_scaffolded -d /home/ubuntu/tools/L_RNA_scaffolder/ -i ${POLISHED%.fasta}.BLAT_ONT.psl -j $POLISHED \
> ${POLISHED%.fasta}_scaffolded/${POLISHED%.fasta}_L_RNA_scaffolding.log 2> ${POLISHED%.fasta}_scaffolded/${POLISHED%.fasta}_L_RNA_scaffolding.err ; echo "done" ; done
