#!/bin/bash
#
#SBATCH --job-name=LRNA_GMAP
#SBATCH --cluster=cm2_tiny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user m.eitel@lmu.de
#SBATCH --mail-type=END,FAIL



NUMTH=28
READS=/dss/dssfs02/lwp-dss-0001/uk213/uk213-dss-0000/di29vos/avas/data/RNA/ONT/Aphrocallistes_ONT_cDNA_run1-3_lordec_100bp.fasta

# to generate the index of the genome gmap_build
for ASSEMBLY in *hardmasked.fasta ; do gmap_build -d ${ASSEMBLY%.fasta} $ASSEMBLY ; done


# to generate sam file (mapping of the sdna genome)
for ASSEMBLY in *hardmasked.fasta ; do \
gmap -d ${ASSEMBLY%.fasta} -k 15 -B 4 -t $NUMTH \
--cross-species -A --exons=cdna --format=samse --npaths=0 --sam-extended-cigar --nofails --no-chimeras --max-intronlength-middle=50000 --max-intronlength-ends=50000 \
$READS > ${ASSEMBLY%.fasta}.ONT-RNA_gmap.sam ; done


#convert sam to bam and sort
for ASSEMBLY in *hardmasked.fasta ; do samtools view -F 4 -@ $NUMTH -b -S ${ASSEMBLY%.fasta}.ONT-RNA_gmap.sam | \
samtools sort -@ $NUMTH > ${ASSEMBLY%.fasta}.ONT-RNA_gmap.sorted.bam ; done

rm *.sam


# count reads
mkdir counts
for ASSEMBLY in *hardmasked.fasta ; do samtools view  -@ NUMTH -F 0x904 -c ${ASSEMBLY%.fasta}.ONT-RNA_gmap.sorted.bam > counts/${ASSEMBLY%.fasta}.ONT-RNA_gmap.sorted.counts ; done

