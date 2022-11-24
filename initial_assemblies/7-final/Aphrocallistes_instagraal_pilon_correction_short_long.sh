#!/bin/bash
#
#SBATCH --job-name=PILON
#SBATCH --cluster=cm2_tiny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user m.eitel@lmu.de
#SBATCH --mail-type=END,FAIL

NUMTH=28
FRAGS=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.PE-DNA_bowtie2.sorted.bam
NANOPORE=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.ONT-DNA_minimap2.sorted.bam
REF=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted.fasta


samtools index -@ $NUMTH Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.PE-DNA_bowtie2.sorted.bam
samtools index -@ $NUMTH Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.ONT-DNA_minimap2.sorted.bam

pilon --genome $REF \
--frags $FRAGS \
--nanopore $NANOPORE \
--threads $NUMTH \
--vcf \
--changes \
--tracks \
--diploid \
--fix "indels" \
--output ${REF%.fasta}.pilon_short_long \
--outdir ${REF%.fasta}.pilon_short_long \
> Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted.pilon_short_long.log


sed -i 's/_pilon//' *fasta