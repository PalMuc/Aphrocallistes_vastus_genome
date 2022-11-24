#!/bin/bash


R1=/home/ubuntu/avas/data/DNA/PE/karect_Avastus_S0_L001_R1.t.f.fastq.gz
R2=/home/ubuntu/avas/data/DNA/PE/karect_Avastus_S0_L001_R2.t.f.fastq.gz
LONGR=/home/ubuntu/avas/data/DNA/ONT/porechop/Aphrocallistes_20171103_1624_MinION_DNA_guppy-3.4.5_R9.4.1_450bps_hac_porechopped_1kb.fasta.gz
LRTYPE=map-ont
NUMTH=40
GSIZE=80m


# mapping corrected (karect) Illumina PE reads
#for DRAFT in *fasta ; do echo "minimap2 --secondary=no --MD -ax sr -t $NUMTH $DRAFT $R1 $R2 | samtools view -Sb - > ${DRAFT%.fasta}.mapped_PE.bam" ; \
#minimap2 --secondary=no --MD -ax sr -t $NUMTH $DRAFT $R1 $R2 2>> minimap2_PE.log | samtools view -Sb - > ${DRAFT%.fasta}.mapped_PE.bam ; echo "done" ; done

for PEMAP in *.mapped_PE.bam ; do echo "samtools sort -@ $NUMTH $PEMAP > ${PEMAP%.bam}.sorted.bam" ; \
samtools sort -@ $NUMTH $PEMAP > ${PEMAP%.bam}.sorted.bam ; echo "done" ; done

for PEMAPSORT in *.mapped_PE.bam ; do samtools index $PEMAPSORT -@ $NUMTH ; done

#rm *.mapped_PE.bam

# mapping raw ONR reads
for DRAFT in *fasta ; do echo "minimap2 --secondary=no --MD -ax $LRTYPE -t $NUMTH $DRAFT $LONGR | samtools view -Sb - > ${DRAFT%.fasta}.mapped_ONT.bam" ; \
minimap2 --secondary=no --MD -ax $LRTYPE -t $NUMTH $DRAFT $LONGR 2>> minimap2_ONT.log | samtools view -Sb - > ${DRAFT%.fasta}.mapped_ONT.bam ; echo "done" ; done

for ONTMAP in *.mapped_ONT.bam ; do echo "samtools sort -@ $NUMTH $ONTMAP > ${ONTMAP%.bam}.sorted.bam" ; \
samtools sort -@ $NUMTH $ONTMAP > ${ONTMAP%.bam}.sorted.bam ; echo "done" ; done

#rm *.mapped_ONT.bam

# generate indices of mapped reads:
for MAPPED in *.sorted.bam ; do samtools index $MAPPED -@ $NUMTH ; done

# mapping statistics:
mkdir flagstat
for FILE in *.sorted.bam ; do "samtools flagstat -@ $NUMTH $FILE > flagstat/${FILE%.bam}.flagstat" ; \
samtools flagstat -@ $NUMTH $FILE > flagstat/${FILE%.bam}.flagstat ; echo "done" ; done


## running HyPo:
echo -e "$R1\n$R2" > short_reads_names.txt

mkdir output

for DRAFT in *fasta ; do echo "hypo -d $DRAFT -r @short_reads_names.txt -s $GSIZE -c 70 -p 96 -t $NUMTH \
-b ${DRAFT%.fasta}.mapped_PE.sorted.bam -B ${DRAFT%.fasta}.mapped_ONT.sorted.bam -o output/${DRAFT%.fasta}.hypo.fasta" ; \
hypo -d $DRAFT -r @short_reads_names.txt -s $GSIZE -c 70 -p 96 -t $NUMTH \
-b ${DRAFT%.fasta}.mapped_PE.sorted.bam -B ${DRAFT%.fasta}.mapped_ONT.sorted.bam -o output/${DRAFT%.fasta}.hypo.fasta ; echo "done" ; done














