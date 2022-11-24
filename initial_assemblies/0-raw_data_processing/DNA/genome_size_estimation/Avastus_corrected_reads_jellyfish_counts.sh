#!/bin/bash
#
#SBATCH --job-name=GS_AVAS
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mail-user m.eitel@lmu.de
#SBATCH --mail-type=END,FAIL
#SBATCH --err=GS_AVAS.err
#SBATCH --mem=100G




# kmer=21
#zcat /home/meitel/data/avas/data/DNA/PE/karect_Avastus_S0_L001_combined.t.f.fastq.gz | jellyfish count /dev/fd/0 -C -m 21 -s 20G -t 20 -o Avastus_corrected_reads_jellyfish_counts_k21.jf
#jellyfish histo -t 10 Avastus_corrected_reads_jellyfish_counts_k21.jf > Avastus_corrected_reads_jellyfish_counts_k21.histo


# kmer=25
zcat /home/meitel/data/avas/data/DNA/PE/karect_Avastus_S0_L001_combined.t.f.fastq.gz | jellyfish count /dev/fd/0 -C -m 25 -s 10G -t 20 -o Avastus_corrected_reads_jellyfish_counts_k25.jf
jellyfish histo -t 10 Avastus_corrected_reads_jellyfish_counts_k25.jf > Avastus_corrected_reads_jellyfish_counts_k25.histo


# kmer=31
zcat /home/meitel/data/avas/data/DNA/PE/karect_Avastus_S0_L001_combined.t.f.fastq.gz | jellyfish count /dev/fd/0 -C -m 31 -s 10G -t 20 -o Avastus_corrected_reads_jellyfish_counts_k31.jf
jellyfish histo -t 10 Avastus_corrected_reads_jellyfish_counts_k31.jf > Avastus_corrected_reads_jellyfish_counts_k31.histo

