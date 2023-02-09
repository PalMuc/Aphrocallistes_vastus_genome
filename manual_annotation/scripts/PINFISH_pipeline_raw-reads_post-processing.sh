#!/bin/bash


# run on lorean1 VM

###############
### PRESETS ###
###############



NUMTH=20
REF=Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.fasta



###################################
### TRANSDECODER $ REFormatting ###
###################################


gffread -T PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed.gff > PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed.gtf

# extract STRINGTIE transcripts and generate a gff3 file:
gtf_genome_to_cdna_fasta.pl PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed.gtf $REF > PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed_mRNA.fasta
gtf_to_alignment_gff3.pl PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed.gtf > PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed.gff3

# predict coding sequences with TRANSDECODER v.2.1.0:
TransDecoder.LongOrfs -t PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed_mRNA.fasta
TransDecoder.Predict -t PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed_mRNA.fasta


#mv and rename transdecoder outputs:
mv PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed_mRNA.fasta.transdecoder.cds PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed_transdecoder_CDS.fasta
mv PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed_mRNA.fasta.transdecoder.gff3 PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed_transdecoder_mRNA.gff3
mv PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed_mRNA.fasta.transdecoder.pep PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed_transdecoder_proteins.fasta


# generate STRINGTIE transdecoder CDS .gff3 file with GMAP:
#gmap_build -t $NUMTH -d ${REF%.fasta} $REF
gmap -d ${REF%.fasta} PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed_transdecoder_CDS.fasta -f 3 -B 4 -t $NUMTH -n 1 > PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed_transdecoder_CDS.gff3 2> PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed_transdecoder_CDS.log

# generate STRINGTIE transdecoder mRNA .gff3 file with GMAP:
gmap -d ${REF%.fasta} PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed_mRNA.fasta -f 3 -B 4 -t $NUMTH -n 1 > PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed_mRNA.gff3 2> PINFISH_pipeline_raw-reads_clustered_transcripts_collapsed_mRNA.log


conda deactivate


