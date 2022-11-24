#!/bin/sh

conda activate genomics

# Aphrocallistes  cDNA GPU-Guppy-3.4.4 basecalling:

full_dir="/home/ubuntu/avas/data/RNA/ONT/raw/run2_R10/fast5_fail"
taxon="Aphrocallistes"
version="3.4.4"
chemistry="R10"
guppy_basecaller -r -x cuda:0 -i $full_dir -c dna_r10_450bps_hac.cfg --compress_fastq --save_path ${taxon}-cDNA-guppy-${version}-${chemistry}-450bps_hac_fast5_fail

conda deactivate