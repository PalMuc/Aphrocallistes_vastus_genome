#!/bin/bash

NUMTH=40

P1R1=Avastus_PE150_R1.fastq
P1R2=Avastus_PE150_R2.fastq
P2R1=Avastus_PE250_R1.fastq
P2R2=Avastus_PE250_R2.fastq

karect -correct \
-threads=$NUMTH \
-celltype=diploid \
-matchtype=hamming \
-inputfile=$P1R1 \
-inputfile=$P1R2 \
-inputfile=$P2R1 \
-inputfile=$P2R2


