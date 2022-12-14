#!/bin/bash


REF=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.fasta

BAM1=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.PE-RNA_hisat2.sorted.bam
BAM2=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.PE-RNA_hisat2_stranded.fwd.bam
BAM3=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.PE-RNA_hisat2_stranded.rev.bam
BAM4=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.sorted.bam


#sizecutter.py -f $REF | sed s/" "/"	"/ > ${REF%.fasta}.sizes

# PE-RNA_hisat2
bedtools genomecov -ibam $BAM1 -bg -split > ${BAM1%.bam}.bg
bedSort ${BAM1%.bam}.bg ${BAM1%.bam}.bg
bedGraphToBigWig ${BAM1%.bam}.bg ${REF%.fasta}.sizes ${BAM1%.bam}.bg.cov.bw

# PE-RNA_hisat2_stranded_fwd
bedtools genomecov -ibam $BAM2 -bg -split > ${BAM2%.bam}.bg
bedSort ${BAM2%.bam}.bg ${BAM2%.bam}.bg
bedGraphToBigWig ${BAM2%.bam}.bg ${REF%.fasta}.sizes ${BAM2%.bam}.bg.cov.bw

# PE-RNA_hisat2_stranded_rev
bedtools genomecov -ibam $BAM3 -bg -split > ${BAM3%.bam}.bg
bedSort ${BAM3%.bam}.bg ${BAM3%.bam}.bg
bedGraphToBigWig ${BAM3%.bam}.bg ${REF%.fasta}.sizes ${BAM3%.bam}.bg.cov.bw

# ONT-RNA_minimap2
bedtools genomecov -ibam $BAM4 -bg -split > ${BAM4%.bam}.bg
bedSort ${BAM4%.bam}.bg ${BAM4%.bam}.bg
bedGraphToBigWig ${BAM4%.bam}.bg ${REF%.fasta}.sizes ${BAM4%.bam}.bg.cov.bw
