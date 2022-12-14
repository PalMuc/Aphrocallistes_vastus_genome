#!/bin/bash


REF=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.fasta



BAM1=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.PE-DNA_hisat2.bam
BAM2=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.PE-DNA_bowtie2.sorted.bam
BAM3=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.ONT-DNA_minimap2.sorted.bam


#sizecutter.py -f $REF | sed s/" "/"	"/ > ${REF%.fasta}.sizes

# PE-DNA_hisat2
bedtools genomecov -ibam $BAM1 -bg -split > ${BAM1%.bam}.bg
bedSort ${BAM1%.bam}.bg ${BAM1%.bam}.bg
bedGraphToBigWig ${BAM1%.bam}.bg ${REF%.fasta}.sizes ${BAM1%.bam}.bg.cov.bw

# PE-DNA_bowtie2
bedtools genomecov -ibam $BAM2 -bg -split > ${BAM2%.bam}.bg
bedSort ${BAM2%.bam}.bg ${BAM2%.bam}.bg
bedGraphToBigWig ${BAM2%.bam}.bg ${REF%.fasta}.sizes ${BAM2%.bam}.bg.cov.bw

# ONT-DNA_minimap2
bedtools genomecov -ibam $BAM3 -bg -split > ${BAM3%.bam}.bg
bedSort ${BAM3%.bam}.bg ${BAM3%.bam}.bg
bedGraphToBigWig ${BAM3%.bam}.bg ${REF%.fasta}.sizes ${BAM3%.bam}.bg.cov.bw


mv *.bg bedgraph
mv *.bw bigwigs