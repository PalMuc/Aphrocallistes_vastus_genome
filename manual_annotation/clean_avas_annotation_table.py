#!/usr/bin/env python
#
# clean_avas_annotation_table.py

'''clean_avas_annotation_table.py  last modified 2021-06-22
    for A vastus manual annotation
    removes some unused columns, cleans up formats

./clean_avas_annotation_table.py Avas_manual_annotation_2021-06-22.tsv > Avas_manual_annotation_2021-06-22.fix.tsv

   then make annotation with:

~/git/genome-reannotations/remake_genome_annotation.py -m manually_edited_genes.gtf -i Avas_manual_annotation_2021-06-24.fix.tsv -g BRAKER2_ONT-RNA_minimap2_augustus.unique.w_flags.gff3 -a STRG_ONT-RNA_minimap2_f-0.3.unique.gff PINFISH_pipeline_raw-reads_clustered_transcripts.unique.gff BRAKER2_ref_protein_alignment_gth_w_scafID.unique.gff3 BRAKER2_PE-RNA_augustus.hints.unique.gff3 STRG_PE-RNA_hisat2_stranded_f-0.3.unique.gff STRG_PE-RNA_hisat2_non-stranded_f-0.3.unique.gff PINFISH_stepwise_5prime-TSL-trimmed_stranded_clustered_transcripts.unique.gff PINFISH_stepwise_corrected-reads_clustered_transcripts.unique.gff -s Avas.v1.a0.9 -f Avas.v1.a0.9_names.tab > Avas.v1.a0.9_annotations.gff

'''

import sys
from collections import defaultdict

if len(sys.argv) < 2:
	sys.exit(__doc__)
else:
	raw_annot_table = sys.argv[1]
	frameshift_cds_filename = "{}.fr_cds.fa".format( raw_annot_table.rsplit(".",1)[1] )
	sys.stderr.write("# removing comment lines, header line, was_checked column, and any columns after comment_for_annotators\n")
	linecounter = 0
	annot_counts = defaultdict(int)
	for line in open(sys.argv[1],'r'):
		if line.strip() and line[0]!="#":
			lsplits = line.split("\t")
			if lsplits[0]=="full_scaf_name":
				continue

			annotator = lsplits[2]
			annot_counts[annotator] += 1

			# rename scaffolds
			base_scaffold = lsplits[0]
			if base_scaffold.strip():
				full_scaf_name = "Aphrocalllistes_vastus_HiC-{}".format(base_scaffold)
				lsplits[0] = full_scaf_name
			else: # reuse most recent name
				lsplits[0] = full_scaf_name

			
			if lsplits[9].find("CDS")>-1:
				if lsplits[8].strip(): #
					lsplits[8] = "{}, CDS manually fixed due to error in assembly".format( lsplits[8] )
				else:
					lsplits[8] = "CDS manually fixed due to error in assembly"

			linecounter += 1
			newsplits = lsplits[0:2] + lsplits[3:8] + [""] + lsplits[8:9]
			outline = "{}\n".format( "\t".join(newsplits) )
			sys.stdout.write( outline )

	sys.stderr.write("# counted {} lines\n".format(linecounter) )
	sys.stderr.write("# annotations by:\n".format(linecounter) )
	for k,v in sorted(annot_counts.items(), key=lambda x: x[1], reverse=True):
		sys.stderr.write("# {} : {}\n".format(k, v) )
