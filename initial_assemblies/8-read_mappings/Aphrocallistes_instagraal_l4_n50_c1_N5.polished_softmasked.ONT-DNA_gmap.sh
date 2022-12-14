#!/bin/bash
#
#SBATCH --job-name=LDNA_GMAP
#SBATCH --cluster=cm2_tiny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user m.eitel@lmu.de
#SBATCH --mail-type=END,FAIL





NUMTH=28
READS=/dss/dssfs02/lwp-dss-0001/uk213/uk213-dss-0000/di29vos/avas/data/DNA/ONT/Aphrocallistes_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_lordec_combined_1kb.fasta



# to generate the index of the genome gmap_build
for ASSEMBLY in *softmasked.fasta ; do gmap_build -d ${ASSEMBLY%.fasta} $ASSEMBLY ; done


# to generate sam file (mapping of the sdna genome)
for ASSEMBLY in *softmasked.fasta ; do gmap -d ${ASSEMBLY%.fasta} -k 15 -A -B 4 -t $NUMTH \
--cross-species -A --nosplicing --format=samse --npaths=0 --sam-extended-cigar --no-chimeras --nofails $READS \
> ${ASSEMBLY%.fasta}.ONT-DNA_gmap.sam ; done


#convert sam to bam and sort
for ASSEMBLY in *softmasked.fasta ; do samtools view -F 4 -@ 12 -b -S ${ASSEMBLY%.fasta}.ONT-DNA_gmap.sam | samtools sort -@ 12> ${ASSEMBLY%.fasta}.ONT-DNA_gmap.sorted.bam ; done


# count reads
mkdir counts
for ASSEMBLY in *softmasked.fasta ; do samtools view  -@ NUMTH -F 0x904 -c ${ASSEMBLY%.fasta}.ONT-DNA_gmap.sorted.bam > counts/${ASSEMBLY%.fasta}.ONT-DNA_gmap.sorted.counts ; done

