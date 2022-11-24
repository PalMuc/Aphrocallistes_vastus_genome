#!/bin/bash



# LORDEC error correction of 1kb reads with k19:
lordec-correct -T 20 -k 19 -s 3 \
-i /home/ubuntu/avas/data/DNA/ONT/porechop/Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_1kb.fasta.gz \
-2 /home/ubuntu/avas/data/DNA/PE/karect_Avastus_S0_L001_combined.t.f.fastq.gz \
-o Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_1kb_lordec_k19.fasta

gzip Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_1kb_lordec_k19.fasta
rm /home/ubuntu/avas/data/DNA/PE/*h5



# LORDEC error correction of 1kb reads with k33:
lordec-correct -T 20 -k 33 -s 3 \
-i Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_1kb_lordec_k19.fasta.gz  \
-2 /home/ubuntu/avas/data/DNA/PE/karect_Avastus_S0_L001_combined.t.f.fastq.gz \
-o Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_1kb_lordec_k33.fasta

gzip Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_1kb_lordec_k33.fasta
rm /home/ubuntu/avas/data/DNA/PE/*h5




# LORDEC error correction of 1kb reads with k41:
lordec-correct -T 20 -k 41 -s 3 \
-i Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_1kb_lordec_k33.fasta.gz  \
-2 /home/ubuntu/avas/data/DNA/PE/karect_Avastus_S0_L001_combined.t.f.fastq.gz \
-o Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_1kb_lordec_k41.fasta \

rm /home/ubuntu/avas/data/DNA/PE/*h5



# trimming of error corrected of 1kb reads (after 3 iterations with k19,k33,k41):
lordec-trim -i Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_1kb_lordec_k41.fasta \
-o Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_1kb_lordec_k41_trimmed.fasta

gzip Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_1kb_lordec_k41.fasta
gzip Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_1kb_lordec_k41_trimmed.fasta





# size filtering of corrected & trimmed reads:
gunzip -c Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_1kb_lordec_k41_trimmed.fasta.gz | \
sizecutter.py -n -a 1000 -p - 2> Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_1kb.fasta.stats | \
gzip > Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_1kb.fasta.gz

gunzip -c Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_1kb.fasta.gz | \
sizecutter.py -n -a 5000 -p - 2> Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_5kb.fasta.stats | \
gzip > Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_5kb.fasta.gz

gunzip -c Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_5kb.fasta.gz | \
sizecutter.py -n -a 10000 -p - 2> Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_10kb.fasta.stats | \
gzip > Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_10kb.fasta.gz

gunzip -c Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_10kb.fasta.gz | \
sizecutter.py -n -a 15000 -p - 2> Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_15kb.fasta.stats | \
gzip > Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_15kb.fasta.gz

gunzip -c Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_15kb.fasta.gz | \
sizecutter.py -n -a 20000 -p - 2> Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_20kb.fasta.stats | \
gzip > Aphrocallistes_20171109_1116_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_20kb.fasta.gz

