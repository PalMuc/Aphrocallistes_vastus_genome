#!/bin/sh




# Aphrocallistes_20180308_1413_PromethION_DNA_guppy-3.4.5_R9.4.1_450bps_hac.fastq.gz
porechop -t 20 \
-i /home/ubuntu/avas/data/DNA/ONT/raw/Aphrocallistes_20180308_1413_PromethION_DNA_guppy-3.4.5_R9.4.1_450bps_hac.fastq.gz \
-o /home/ubuntu/avas/data/DNA/ONT/porechop/Aphrocallistes_20180308_1413_PromethION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped.fastq.gz \
--adapter_threshold 80.0 \
--check_reads 100000 \
--extra_end_trim 10 \
--end_threshold 70 \
--format fastq.gz

gunzip -c /home/ubuntu/avas/data/DNA/ONT/porechop/Aphrocallistes_20180308_1413_PromethION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped.fastq.gz | \
sizecutter.py -n -i fastq -o fasta -a 1000 -p - \
2> /home/ubuntu/avas/data/DNA/ONT/porechop/Aphrocallistes_20180308_1413_PromethION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_1kb.stats | \
gzip > /home/ubuntu/avas/data/DNA/ONT/porechop/Aphrocallistes_20180308_1413_PromethION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_1kb.fasta.gz

