#!/bin/bash
#
#SBATCH --job-name=GMAP
#SBATCH --cluster=cm2_tiny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user m.eitel@lmu.de
#SBATCH --mail-type=END,FAIL


NUMTH=28
TRANSCRIPTS=Aphrocallistes_RNAseq_Ind_1-5_combined_f25_reheader.combined.okay.fa

# to generate the index of the genome gmap_build
#for ASSEMBLY in *hardmasked.fasta ; do gmap_build -d ${ASSEMBLY%.fasta} $ASSEMBLY ; done


# to generate sam file (mapping of the sdna genome)
for ASSEMBLY in *hardmasked.fasta ; do \
gmap -d ${ASSEMBLY%.fasta} -B 4 -t $NUMTH -O -A --cross-species --no-chimeras --exons=cdna --format=samse --npaths=1 --sam-extended-cigar --nofails \
$TRANSCRIPTS > ${ASSEMBLY%.fasta}.transcriptome_gmap.sam ; done

# convert sam to bam and sort
for ASSEMBLY in *hardmasked.fasta ; do samtools view -F 4 -@ $NUMTH -b -S ${ASSEMBLY%.fasta}.transcriptome_gmap.sam | samtools sort -@ $NUMTH > ${ASSEMBLY%.fasta}.transcriptome_gmap.sorted.bam ; done

# count reads
mkdir counts
for ASSEMBLY in *hardmasked.fasta ; do samtools view  -@ NUMTH -F 0x904 -c ${ASSEMBLY%.fasta}.transcriptome_gmap.sorted.bam > counts/${ASSEMBLY%.fasta}.transcriptome_gmap.sorted.counts ; done

