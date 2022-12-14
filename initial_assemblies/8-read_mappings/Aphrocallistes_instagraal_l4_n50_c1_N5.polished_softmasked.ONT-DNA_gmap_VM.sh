#!/bin/bash



NUMTH=40
READS=/home/ubuntu/avas/data/DNA/ONT/lordec/Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_combined_1kb.fasta.gz
REF=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted.fasta


# to generate the index of the genome gmap_build
gmap_build -d ${REF%.fasta} $REF


# to generate sam file (mapping of the sdna genome)
gmap -d ${REF%.fasta} -k 15 -A -B 4 -t $NUMTH \
--cross-species -A --nosplicing --format=samse --npaths=0 --sam-extended-cigar --no-chimeras --nofails $READS \
> ${REF%.fasta}.ONT-DNA_gmap.sam


#convert sam to bam and sort
samtools view -F 4 -@ $NUMTH -b -S ${REF%.fasta}.ONT-DNA_gmap.sam | samtools sort -@ $NUMTH > ${REF%.fasta}.ONT-DNA_gmap.sorted.bam


# count reads
mkdir counts
samtools view  -@ $NUMTH -F 0x904 -c ${REF%.fasta}.ONT-DNA_gmap.sorted.bam > counts/${REF%.fasta}.ONT-DNA_gmap.sorted.counts

