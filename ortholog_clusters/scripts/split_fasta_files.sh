#!/bin/bash



NUMFILES=2710

mkdir setA ; shuf -n $NUMFILES -e *fasta | xargs -i mv {} setA
mkdir setB ; shuf -n $NUMFILES -e *fasta | xargs -i mv {} setB
mkdir setC ; shuf -n $NUMFILES -e *fasta | xargs -i mv {} setC
mkdir setD ; shuf -n $NUMFILES -e *fasta | xargs -i mv {} setD
mkdir setE ; shuf -n $NUMFILES -e *fasta | xargs -i mv {} setE
mkdir setF ; shuf -n $NUMFILES -e *fasta | xargs -i mv {} setF
mkdir setG ; shuf -n $NUMFILES -e *fasta | xargs -i mv {} setG
mkdir setH ; shuf -n $NUMFILES -e *fasta | xargs -i mv {} setH

