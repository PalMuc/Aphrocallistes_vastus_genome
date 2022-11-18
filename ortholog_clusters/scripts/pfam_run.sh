#!/bin/bash


NUMTH=4
PFAM=/home/ubuntu/avas/ortholog_assignments/makehomologs/pfam_hmm/Pfam-A.hmm

# search PfamA with HMMSEARCH:
for file in *fasta ; do hmmscan --cpu $NUMTH --domtblout ${file%.fasta}.pfam.txt $PFAM $file ; done

# reformat output:
for file in *txt ; do tabularpfam.py $file > ${file%.txt}.tab ; done

mkdir original_output
mv *txt original_output
mkdir formatted_output
mv *tab formatted_output

