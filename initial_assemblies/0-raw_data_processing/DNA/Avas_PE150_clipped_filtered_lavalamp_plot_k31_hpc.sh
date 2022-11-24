#!/bin/bash
#
#SBATCH --job-name=lavAvas
#SBATCH --ntasks=40
#SBATCH --mail-user m.eitel@lmu.de
#SBATCH --mail-type=ALL
#SBATCH --err=Avas_PE150_clipped_filtered_lavalamp_plot_k31.err
#SBATCH --mem=200000



#gunzip -c /home/meitel/data/avas/data/gDNA/Avastus_S0_L001_R1.t.f.fastq.gz > /home/meitel/data/avas/data/gDNA/Avastus_S0_L001_R1.t.f.fastq
#gunzip -c /home/meitel/data/avas/data/gDNA/Avastus_S0_L001_R2.t.f.fastq.gz > /home/meitel/data/avas/data/gDNA/Avastus_S0_L001_R2.t.f.fastq

# step 1:
#cat /home/meitel/data/avas/data/gDNA/Avastus_S0_L001_R1.t.f.fastq /home/meitel/data/avas/data/gDNA/Avastus_S0_L001_R2.t.f.fastq | jellyfish count -m 31 -s 4G -C -U 500 -t 40 -o avas_PE150_clipped_filtered_k31.counts /dev/fd/0

# step 2:
#jellyfish dump avas_PE150_clipped_filtered_k31.counts | fastqdumps2histo.py -k 31 -u 500 -j avas_PE150_clipped_filtered_k31.dumps - > avas_PE150_clipped_filtered_k31_counts_gc_cov_histo.csv

# step 3:
kmersorter.py -T -k 31 -p 40 -D /home/meitel/tools/trinityrnaseq-Trinity-v2.5.1/ -1 /home/meitel/data/avas/data/gDNA/Avastus_S0_L001_R1.t.f.fastq -2 /home/meitel/data/avas/data/gDNA/Avastus_S0_L001_R2.t.f.fastq avas_PE150_clipped_filtered_k31.dumps

# step 4:
fastqdumps2histo.py -s /home/meitel/data/avas/data/gDNA/Avastus_S0_L001_R1.t.f.fasta.k31.stats /home/meitel/data/avas/data/gDNA/Avastus_S0_L001_R2.t.f.fasta.k31.stats -f /home/meitel/data/avas/data/gDNA/Avastus_S0_L001_R1.t.f.fastq /home/meitel/data/avas/data/gDNA/Avastus_S0_L001_R2.t.f.fastq -u 500 -T - > avas_PE150_clipped_filtered_k31_reads_gc_cov_histo.csv