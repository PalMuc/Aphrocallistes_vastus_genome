#!/bin/sh

# Aphrocallistes_20180308_1413_PromethION DNA GPU-Guppy-3.4.5 basecalling:

full_dir="/home/ubuntu/avas/data/DNA/ONT/raw/Aphrocallistes_20180308_1413_PromethION/fast5"
taxon="Aphrocallistes_20180308_1413_PromethION"
version="3.4.5"
chemistry="R9.4.1"



guppy_basecaller -r -x cuda:0,1 -i $full_dir \
-c dna_r9.4.1_450bps_hac_prom.cfg  \
--compress_fastq \
--chunk_size 2000 \
--num_callers 12 \
--chunks_per_runner 4096 \
--save_path ${taxon}_DNA_guppy-${version}_${chemistry}_450bps_hac

cat /home/ubuntu/avas/data/DNA/ONT/raw/Aphrocallistes_20180308_1413_PromethION_DNA_guppy-3.4.5_R9.4.1_450bps_hac/*gz \
> Aphrocallistes_20180308_1413_PromethION_DNA_guppy-3.4.5_R9.4.1_450bps_hac.fastq.gz