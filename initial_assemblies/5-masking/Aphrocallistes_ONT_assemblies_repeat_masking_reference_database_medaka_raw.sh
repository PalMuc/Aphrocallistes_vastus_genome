#!/bin/bash
#
#SBATCH --job-name=R-MOD
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

MASKDIR=/dss/dssfs02/lwp-dss-0001/uk213/uk213-dss-0000/di29vos/avas/repeatmasking/MEDAKA/purged/raw_reads

REFASS=${MASKDIR}/Aphrocallistes_CANU_assembly-1_5kb_raw_reads_Racon-ONT_final_Medaka-ONT_Racon-PE_RNA_scaffolded_final.purged.fasta


# build genomic source reference database:
BuildDatabase -engine ncbi -name aphrocallistes_vastus_medaka_raw ${REFASS}


# generate taxon specific repeat library:
RepeatModeler -database aphrocallistes_vastus_medaka_raw -engine ncbi -pa $NUMTH


cp RM_*/consensi.fa aphrocallistes_vastus-families_medaka_raw.fa
cp RM_*/families.stk aphrocallistes_vastus-families_medaka_raw.stk


