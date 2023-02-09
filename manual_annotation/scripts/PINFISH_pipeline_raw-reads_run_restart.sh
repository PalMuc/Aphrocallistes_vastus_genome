#!/bin/bash
#
#SBATCH --job-name=PINFISH
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --mem=40G
#SBATCH --qos=low
#SBATCH --mail-user m.eitel@lmu.de
#SBATCH --mail-type=END,FAIL



/home/meitel/tools/pipeline-pinfish-analysis/pinfish/polish_clusters/polish_clusters -t 40 \
-a /home/meitel/data/avas/pinfish/scaffolded_HiC/pipeline/output/results/cluster_memberships.tsv \
-c 5 -o /home/meitel/data/avas/pinfish/scaffolded_HiC/pipeline/output/results/polished_transcripts.fas \
/home/meitel/data/avas/pinfish/scaffolded_HiC/pipeline/output/alignments/reads_aln_sorted.bam
