#!/bin/bash
#
#SBATCH --job-name=MINIMAP2
#SBATCH --cluster=cm2_tiny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user m.eitel@lmu.de
#SBATCH --mail-type=END,FAIL




NUMTH=28
TRANSCRIPTS=Aphrocallistes_RNAseq_Ind_1-5_combined_f25_reheader.combined.okay.fa

# to generate sam file (mapping of the dna genome)
for ASSEMBLY in *hardmasked.fasta ; do minimap2 -ax splice -ub--secondary=no -t $NUMTH $ASSEMBLY $TRANSCRIPTS > ${ASSEMBLY%.fasta}.transcriptome_minimap2.sam ; done


# convert sam to bam and sort
for ASSEMBLY in *hardmasked.fasta ; do samtools view -F 4 -@ $NUMTH -b -S ${ASSEMBLY%.fasta}.transcriptome_minimap2.sam | samtools sort -@ $NUMTH > ${ASSEMBLY%.fasta}.transcriptome_minimap2.sorted.bam ; done

# count reads
mkdir counts
for ASSEMBLY in *hardmasked.fasta ; do samtools view  -@ NUMTH -F 0x904 -c ${ASSEMBLY%.fasta}.transcriptome_minimap2.sorted.bam > counts/${ASSEMBLY%.fasta}.transcriptome_minimap2.sorted.counts ; done

