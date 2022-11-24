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

REF=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted.fasta


samtools index -@ $NUMTH Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.PE-DNA_bowtie2.sorted.bam


pilon --genome $REF \
--frags $FRAGS \
--threads $NUMTH \
--vcf \
--changes \
--tracks \
--diploid \
--fix "indels" \
--output ${REF%.fasta}.pilon_short \
--outdir ${REF%.fasta}.pilon_short \
> Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted.pilon_short.log

sed -i 's/_pilon//' *fasta