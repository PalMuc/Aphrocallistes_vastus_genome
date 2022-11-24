#!/bin/bash
#
#SBATCH --job-name=lavAvas
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user m.eitel@lmu.de
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp2
#SBATCH --time=12:00:00
#SBATCH --mem=50000


module load java/1.8


gunzip Avastus_S0_L001_R1.t.f.fastq.gz
gunzip Avastus_S0_L001_R2.t.f.fastq.gz

# step 1:
cat Avastus_S0_L001_R1.t.f.fastq Avastus_S0_L001_R2.t.f.fastq | /naslx/projects/uk213/di29vos/trinityrnaseq-Trinity-v2.4.0/trinity-plugins/jellyfish/bin/jellyfish count -m 31 -s 4G -C -U 500 -t 28 -o avas_PE150_clipped_filtered_k31.counts /dev/fd/0

# step 2:
/naslx/projects/uk213/di29vos/trinityrnaseq-Trinity-v2.4.0/trinity-plugins/jellyfish/bin/jellyfish dump avas_PE150_clipped_filtered_k31.counts | /naslx/projects/uk213/di29vos/lavaLampPlot-master/fastqdumps2histo.py -k 31 -u 500 -j avas_PE150_clipped_filtered_k31.dumps - > avas_PE150_clipped_filtered_k31_counts_gc_cov_histo.csv

# step 3:
/naslx/projects/uk213/di29vos/lavaLampPlot-master/kmersorter.py -T -k 31 -p 28 -D /naslx/projects/uk213/di29vos/trinityrnaseq-Trinity-v2.4.0/ -1 Avastus_S0_L001_R1.t.f.fastq -2 Avastus_S0_L001_R2.t.f.fastq avas_PE150_clipped_filtered_k31.dumps

# step 4:
/naslx/projects/uk213/di29vos/lavaLampPlot-master/fastqdumps2histo.py -s Avastus_S0_L001_R1.t.f.left.fa.k31.stats Avastus_S0_L001_R2.t.f.right.fa.k31.stats -f Avastus_S0_L001_R1.t.f.fastq Avastus_S0_L001_R2.t.f.fastq -u 500 -T - > avas_PE150_clipped_filtered_k31_reads_gc_cov_histo.csv