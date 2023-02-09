#!/bin/bash

DATA=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.PE-RNA_hisat2_stranded.sorted.bam


# Forward strand.
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand
#
samtools view -@ 8 -b -f 128 -F 16 $DATA > fwd1.bam
samtools index -@ 8 fwd1.bam

samtools view -@ 8 -b -f 80 $DATA > fwd2.bam
samtools index -@ 8 fwd2.bam

#
# Combine alignments that originate on the forward strand.
#
samtools merge -@ 8 -f Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.PE-RNA_hisat2_stranded.fwd.bam fwd1.bam fwd2.bam
samtools index -@ 8 Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.PE-RNA_hisat2_stranded.fwd.bam

# Reverse strand
#
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
#
samtools view -@ 8 -b -f 144 $DATA > rev1.bam
samtools index -@ 8 rev1.bam

samtools view -@ 8 -b -f 64 -F 16 $DATA > rev2.bam
samtools index -@ 8 rev2.bam

#
# Combine alignments that originate on the reverse strand.
#
samtools merge -@ 8 -f Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.PE-RNA_hisat2_stranded.rev.bam rev1.bam rev2.bam
samtools index -@ 8 Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.PE-RNA_hisat2_stranded.rev.bam