#!/bin/bash

NUMTH=4


# create MAFFT-einsi alignments
for file in *fasta ; do mafft --maxiterate 1000 --genafpair --thread $NUMTH $file >$file.aln ; done

# create FASTREE gene trees:
for file in *aln ; do fasttree -wag -gamma $file > $file.tre ; done

