#!/bin/sh

source activate genomics

# Aphrocallistes  cDNA GPU-Guppy-3.2.4 basecalling:

full_dir="/home/ubuntu/avas/data/RNA/ONT/raw/run1_R9.4/Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_fast5_fail"
taxon="Aphrocallistes"
version="3.2.4"
chemistry="R9.4"
guppy_basecaller -r -i $full_dir --save_path ${taxon}-cDNA-guppy-${version}-${chemistry}-450bps_hac_fast5_fail \
-x cuda:1 \
-c dna_r9.4.1_450bps_hac.cfg \
--compress_fastq \
--chunk_size 4000 \
--chunks_per_runner 2048 \
--gpu_runners_per_device 16 \
--num_callers 16


conda deactivate