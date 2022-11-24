#!/bin/bash

source activate genomics


# 1. error correction with k19
lordec-correct -T 20 -k 19 -s 3 \
-i /home/ubuntu/avas/data/RNA/ONT/pychopper/Aphrocallistes-cDNA-guppy-3.4.4-R9.4.1-450bps_hac_run3.pychopper2_all_100bp_to_20kb.fasta.gz \
-2 \
/home/ubuntu/avas/data/RNA/PE/Avastus_RNAseq_LEY_WOE_concat_1.f25.fastq.gz \
/home/ubuntu/avas/data/RNA/PE/Avastus_RNAseq_LEY_WOE_concat_2.f25.fastq.gz \
-o Aphrocallistes-cDNA-guppy-3.4.4-R9.4.1_run3_combined_hac.lordec_k19.fasta \
2> Aphrocallistes-cDNA-guppy-3.4.4-R9.4.1_run3_combined_hac.lordec_k19.log

gzip Aphrocallistes-cDNA-guppy-3.4.4-R9.4.1_run3_combined_hac.lordec_k19.fasta

rm /home/ubuntu/avas/data/RNA/PE/*h5

# 2. error correction with k31
lordec-correct -T 20 -k 31 -s 3 \
-i Aphrocallistes-cDNA-guppy-3.4.4-R9.4.1_run3_combined_hac.lordec_k19.fasta.gz \
-2 \
/home/ubuntu/avas/data/RNA/PE/Avastus_RNAseq_LEY_WOE_concat_1.f25.fastq.gz \
/home/ubuntu/avas/data/RNA/PE/Avastus_RNAseq_LEY_WOE_concat_2.f25.fastq.gz \
-o Aphrocallistes-cDNA-guppy-3.4.4-R9.4.1_run3_combined_hac.lordec_k31.fasta \
2> Aphrocallistes-cDNA-guppy-3.4.4-R9.4.1_run3_combined_hac.lordec_k31.log

gzip Aphrocallistes-cDNA-guppy-3.4.4-R9.4.1_run3_combined_hac.lordec_k31.fasta

rm /home/ubuntu/avas/data/RNA/PE/*h5


# 3. error correction with k41
lordec-correct -T 20 -k 41 -s 3 \
-i Aphrocallistes-cDNA-guppy-3.4.4-R9.4.1_run3_combined_hac.lordec_k31.fasta.gz \
-2 \
/home/ubuntu/avas/data/RNA/PE/Avastus_RNAseq_LEY_WOE_concat_1.f25.fastq.gz \
/home/ubuntu/avas/data/RNA/PE/Avastus_RNAseq_LEY_WOE_concat_2.f25.fastq.gz \
-o Aphrocallistes-cDNA-guppy-3.4.4-R9.4.1_run3_combined_hac.lordec_k41.fasta \
2> Aphrocallistes-cDNA-guppy-3.4.4-R9.4.1_run3_combined_hac.lordec_k41.log

rm /home/ubuntu/avas/data/RNA/PE/*h5


# trimming of uncorrect read ends:
lordec-trim -i Aphrocallistes-cDNA-guppy-3.4.4-R9.4.1_run3_combined_hac.lordec_k41.fasta \
-o Aphrocallistes-cDNA-guppy-3.4.4-R9.4.1_run3_combined_hac.lordec_k41_trimmed.fasta

gzip Aphrocallistes-cDNA-guppy-3.4.4-R9.4.1_run3_combined_hac.lordec_k41.fasta
gzip Aphrocallistes-cDNA-guppy-3.4.4-R9.4.1_run3_combined_hac.lordec_k41_trimmed.fasta

gunzip -c Aphrocallistes-cDNA-guppy-3.4.4-R9.4.1_run3_combined_hac.lordec_k41_trimmed.fasta.gz | \
sizecutter.py -a 100 -p - | gzip > Aphrocallistes-cDNA-guppy-3.4.4-R9.4.1_run3_combined_hac.lordec_100bp.fasta.gz

for file in *fasta.gz; do gunzip -c $file | sizecutter.py -qn - 2> $file.stats ; done
for file in *fasta.gz; do gunzip -c $file | lowercount.pl - > $file.lowercount ; done

conda deactivate


