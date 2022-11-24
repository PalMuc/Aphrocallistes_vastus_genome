#!/bin/bash

###############
### PRESETS ###
###############


## NOTE1: masking taxon-specific repeats following guidlines by "https://blaxter-lab-documentation.readthedocs.io/en/latest/repeat-masking.html"


source activate masking


NUMTH=20

REF=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted.fasta


# build genomic source reference database:
BuildDatabase -engine ncbi -name Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted $REF


# generate taxon specific repeat library:
RepeatModeler -database Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted -engine ncbi -pa $NUMTH

cp RM_*/consensi.fa Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted-families.fa
cp RM_*/families.stk Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted-families.stk

rm -r RM_*

conda deactivate
