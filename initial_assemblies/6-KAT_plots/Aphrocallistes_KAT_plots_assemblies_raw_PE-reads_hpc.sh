#!/bin/bash

NUMTH=12

READS1=Avastus_S0_L001_R1.fastq.gz 
READS2=Avastus_S0_L001_R2.fastq.gz 

# matrix generation:
for file in *fasta; \
do \
kat comp -t $NUMTH -H 1600000000 \
-o ${file%.fasta}.raw_PE-reads \
--output_hists \
--output_type pdf \
'$READS1 $READS2' $file \
> ${file%.fasta}.kat.comp.raw_PE-reads.log ; \
done

# Spectra CN plot:
for file in *fasta; \
do \
kat plot spectra-cn \
-o ${file%.fasta}.kat.spectra-cn_raw_PE-reads \
-x 250 -p pdf \
-t ${file%.fasta}.kat.spectra-cn \
--dpi 600 \
${file%.fasta}.raw_PE-reads-main.mx ; \
done


