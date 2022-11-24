#!/bin/bash


source activate pychopper


cdna_classifier.py -t 10 \
-Y 100000 \
-b /home/ubuntu/tools/anaconda3/envs/pychopper/bin/cDNA_SSP_VNP.fas \
-r Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1_pychopper2_report.pdf \
-u Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2.unclassified.fastq \
-S Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2.stats \
-A Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2.scores \
-x PCS109 \
-B 100000 \
-w Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2.rescued.fastq \
/home/ubuntu/avas/data/RNA/ONT/raw/Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.fastq \
Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2.full_length_output.fastq \
2> Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2.log


conda deactivate

gzip *fastq

source activate genomics

gunzip -c Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2.full_length_output.fastq.gz | \
sizecutter.py -n -i fastq -o fastq -p -a 100 -b 20000 - 2> Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2.full_length_output_100bp_to_20kb.fastq.stats | \
gzip > Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2.full_length_output_100bp_to_20kb.fastq.gz


gunzip -c Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2.rescued.fastq.gz | \
sizecutter.py -n -i fastq -o fastq -p -a 100 -b 20000 - 2> Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2.rescued_100bp_to_20kb.fastq.stats | \
gzip > Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2.rescued_100bp_to_20kb.fastq.gz


gunzip -c Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2.unclassified.fastq.gz | \
sizecutter.py -n -i fastq -o fastq -p -a 100 -b 20000 - 2> Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2.unclassified_100bp_to_20kb.fastq.stats | \
gzip > Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2.unclassified_100bp_to_20kb.fastq.gz

cat *_100bp_to_20kb.fastq.gz > Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2_all_100bp_to_20kb.fastq.gz

gunzip -c Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2_all_100bp_to_20kb.fastq.gz | \
sizecutter.py -n -i fastq -p - 2> Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2_all_100bp_to_20kb.fasta.stats | \
awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' - | \
gzip > Aphrocallistes-cDNA-guppy-3.2.4-R9.4.1-450bps_hac_run1.pychopper2_all_100bp_to_20kb.fasta.gz

conda deactivate





