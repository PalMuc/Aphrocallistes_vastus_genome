#!/bin/bash


# https://github.com/CAFS-bioinformatics/L_RNA_scaffolder


READS=Aphrocallistes-cDNA-guppy-3.4.4-R9.4.1_run3_combined_hac.lordec_100bp.fasta

NUMTH=30

# generate .psl BLAT alignment file:
for POLISHED in *fasta ; do echo "pblat -t=dna -q=dna -noHead -minIdentity=95 -maxIntron=50000 -threads=$NUMTH $POLISHED $CDNAREADS ${POLISHED%.fasta}.BLAT_ONT.psl" ; \
pblat -t=dna -q=dna -noHead -minIdentity=95 -maxIntron=50000 -threads=$NUMTH $POLISHED $READS ${POLISHED%.fasta}.BLAT_ONT.psl ; echo "done" ; done


# running L_RNA_scaffolder
for POLISHED in *.fasta ; do mkdir ${POLISHED%.fasta}_scaffolded ; \
echo "/home/ubuntu/tools/L_RNA_scaffolder/L_RNA_scaffolder.sh -f 2 -o ${POLISHED%.fasta}_scaffolded -d /home/ubuntu/tools/L_RNA_scaffolder/ -i ${POLISHED%.fasta}.BLAT_ONT.psl -j $POLISHED" ; \
/home/ubuntu/tools/L_RNA_scaffolder/L_RNA_scaffolder.sh -f 2 -o ${POLISHED%.fasta}_scaffolded -d /home/ubuntu/tools/L_RNA_scaffolder/ -i ${POLISHED%.fasta}.BLAT_ONT.psl -j $POLISHED \
> ${POLISHED%.fasta}_scaffolded/${POLISHED%.fasta}_L_RNA_scaffolding.log 2> ${POLISHED%.fasta}_scaffolded/${POLISHED%.fasta}_L_RNA_scaffolding.err ; echo "done" ; done
