#!/bin/bash

###############
### PRESETS ###
###############


## NOTE1: masking taxon-specific repeats following guidlines by "https://blaxter-lab-documentation.readthedocs.io/en/latest/repeat-masking.html"
## NOTE2: reference transcriptome is the final TransPi pipeline transdecoder output (TransPi run with kmer setting "C")
## NOTE3: reference proteome is the final TransPi pipeline transdecoder output (TransPi run with kmer setting "C")
## NOTE4: reference repeat database were generated based on REFASS @ Aphrocallistes_ONT_assemblies_repeat_masking_reference_database.sh


source activate masking

mkdir soft
mkdir softmasked_final

mkdir hard
mkdir hardmasked_final

NUMTH=40

REFDB=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted
TRANSCRIPTOME=Aphrocallistes_RNAseq_Ind_1-5_combined_f25_reheader.combined.okay.fa
PROTEOME=Aphrocallistes_RNAseq_Ind_1-5_combined_f25_reheader.combined.okay.fa.transdecoder.pep



# Filter RepeatModeler Library for non-TE protein coding sequences:
## BLAST TRANSPI proteins of all transcriptomes (A-G+S+Z) against RepeatMasker TE database:

blastp -query $PROTEOME \
-db /home/ubuntu/tools/anaconda3/envs/masking/share/RepeatMasker/Libraries/RepeatPeps.lib \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-max_target_seqs 25 \
-culling_limit 2 \
-num_threads $NUMTH \
-evalue 1e-5 \
-out ${PROTEOME}.proteins.fa.vs.RepeatPeps.25cul2.1e5.blastp.out


## Remove TEs from proteome:
fastaqual_select.pl -f $TRANSCRIPTOME -e <(awk '{print $1}' ${PROTEOME}.proteins.fa.vs.RepeatPeps.25cul2.1e5.blastp.out | sort | uniq) \
> ${TRANSCRIPTOME%.fa}_no_tes.fasta


## Blast transcriptome against RepeatModeler library:
makeblastdb -in ${TRANSCRIPTOME%.fa}_no_tes.fasta -dbtype nucl

blastn -task megablast \
-query ${REFDB}.fasta \
-db ${TRANSCRIPTOME%.fa}_no_tes.fasta \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-max_target_seqs 25 \
-culling_limit 2 \
-num_threads $NUMTH \
-evalue 1e-10 \
-out aphrocallistes_vastus_repeatmodeler_lib.vs.transcripts.no_tes.25cul2.1e10.megablast.out


## Remove hits from RepeatModeler library
fastaqual_select.pl -f ${REFDB}.fasta \
-e <(awk '{print $1}' aphrocallistes_vastus_repeatmodeler_lib.vs.transcripts.no_tes.25cul2.1e10.megablast.out | sort | uniq) > \
${REFDB}.filtered_for_CDS_repeats.fa


## repeat masking using combined taxon specific and repeatmodeler library:
# soft masking:
for ASSEMBLY in *fasta ; do \
RepeatMasker -s -pa $NUMTH -norna -dir ./soft -lib ${REFDB}.filtered_for_CDS_repeats.fa \
-a -inv -lcambig -xsmall -gff $ASSEMBLY ; \
done

cp ./soft/*.masked ./softmasked_final/
cd softmasked_final
ls *.masked | awk '{print("mv "$1" "$1)}' | sed 's/\.fasta\.masked/_softmasked\.fasta/2' | /bin/sh
cd ..


# hard masking:
for ASSEMBLY in *fasta ; do \
RepeatMasker -s -pa $NUMTH -norna -dir ./hard -lib ${REFDB}.filtered_for_CDS_repeats.fa \
-a -inv -lcambig -gff $ASSEMBLY ; \
done

cp ./hard/*.masked ./hardmasked_final/
cd hardmasked_final
ls *.masked | awk '{print("mv "$1" "$1)}' | sed 's/\.fasta\.masked/_hardmasked\.fasta/2' | /bin/sh
cd ..


conda deactivate
