#!/bin/bash
#
#SBATCH --job-name=C-MOD
#SBATCH --cluster=cm2_tiny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user m.eitel@lmu.de
#SBATCH --mail-type=END,FAIL


###############
### PRESETS ###
###############


## NOTE1: masking taxon-specific repeats following guidlines by "https://blaxter-lab-documentation.readthedocs.io/en/latest/repeat-masking.html"



NUMTH=28

MASKDIR=/dss/dssfs02/lwp-dss-0001/uk213/uk213-dss-0000/di29vos/avas/repeatmasking/HYPO/purged/corrected_reads

REFASS=${MASKDIR}/Aphrocallistes_CANU_assembly-1_5kb_corrected_reads.hypo_L_RNA_scaffolded_ONT_TR_PE_sorted.purged.fasta


# build genomic source reference database:
BuildDatabase -engine ncbi -name aphrocallistes_vastus_hypo_corrected ${REFASS}


# generate taxon specific repeat library:
RepeatModeler -database aphrocallistes_vastus_hypo_corrected -engine ncbi -pa $NUMTH


cp RM_*/consensi.fa aphrocallistes_vastus-families_hypo_corrected.fa
cp RM_*/families.stk aphrocallistes_vastus-families_hypo_corrected.stk


