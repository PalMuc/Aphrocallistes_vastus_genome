#!/bin/bash
#
#SBATCH --job-name=LRNA_MINIMAP
#SBATCH --cluster=cm2_tiny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user m.eitel@lmu.de
#SBATCH --mail-type=END,FAIL



NUMTH=28
READS=/dss/dssfs02/lwp-dss-0001/uk213/uk213-dss-0000/di29vos/avas/data/RNA/ONT/Aphrocallistes_ONT_cDNA_run1-3_lordec_100bp.fasta

REF=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted.pilon_short_long.fasta

# to generate sam file (mapping of the dna genome)
minimap2 -ax splice -ub -Y --secondary=no -t $NUMTH $REF $READS > ${REF%.fasta}.ONT-RNA_minimap2.sam

# convert sam to bam and sort
samtools view -F 4 -@ 12 -b -S ${REF%.fasta}.ONT-RNA_minimap2.sam | samtools sort -@ 12 > ${REF%.fasta}.ONT-RNA_minimap2.sorted.bam


# count reads
mkdir counts
samtools view  -@ $NUMTH -F 0x904 -c ${REF%.fasta}.ONT-RNA_minimap2.sorted.bam > counts/${REF%.fasta}.ONT-RNA_minimap2.sorted.counts

