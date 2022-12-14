#!/bin/bash
#
#SBATCH --job-name=LRNA_GMAP
#SBATCH --cluster=cm2_tiny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user m.eitel@lmu.de
#SBATCH --mail-type=END,FAIL



NUMTH=28
READS=/dss/dssfs02/lwp-dss-0001/uk213/uk213-dss-0000/di29vos/avas/data/RNA/ONT/Aphrocallistes_ONT_cDNA_run1-3_lordec_TSL_porechopped.100bp.fasta

REF=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted.pilon_short_long.fasta



# to generate the index of the genome gmap_build
gmap_build -d ${REF%.fasta} $REF


# to generate sam file (mapping of the sdna genome)
gmap -d ${REF%.fasta}_k16 -k 16 -B 4 -t $NUMTH -A --exons=cdna --format=samse --npaths=0 --sam-extended-cigar --no-chimeras --split-large-introns  --nofails --max-intronlength-middle=40000 --max-intronlength-ends=40000 \
$READS > ${REF%.fasta}.ONT-RNA_gmap.redone.k16.sam 2> ${REF%.fasta}.ONT-RNA_gmap.redone.log


#convert sam to bam and sort
samtools view -F 4 -@ $NUMTH -b -S ${REF%.fasta}.ONT-RNA_gmap.redone.sam | \
samtools sort -@ $NUMTH > ${REF%.fasta}.ONT-RNA_gmap.sorted.redone.k16.bam
