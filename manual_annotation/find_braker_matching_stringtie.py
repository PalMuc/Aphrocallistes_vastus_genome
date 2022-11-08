#!/usr/bin/env python
#
# find_braker_matching_stringtie.py

'''
find_braker_matching_stringtie.py  last modified 2020-09-15

./find_braker_matching_stringtie.py -s manual_curation_tracks/STRG_ONT-RNA_minimap2_f-0.3.unique.gff -p manual_curation_tracks/PINFISH_pipeline_raw-reads_clustered_transcripts.unique.gff -a manual_curation_tracks/BRAKER2_ONT-RNA_minimap2_augustus.unique.gff3 > manual_curation_tracks/BRAKER2_ONT-RNA_minimap2_augustus.unique.w_flags.gff3


'''

import sys
import time
import argparse
import gzip
from collections import defaultdict

def parse_braker(brakerfile, lr_scaf_tx_dict, lr_tx_exon_dict):
	'''parse braker, check each gene'''

	if brakerfile.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Reading GFF file {} as gzipped\n".format(brakerfile) )
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Reading GFF file {}\n".format(brakerfile) )

	gene_buffer = defaultdict(list) # key is gene id, value is strings of gff lines
	gene_order = []
#HiC-scaffold_001	AUGUSTUS	gene	2844875	2846140	.	+	.	ID=brakerONT.jg3962;Name=brakerONT.jg3962
#HiC-scaffold_001	AUGUSTUS	mRNA	2844875	2846140	.	+	.	ID=brakerONT.jg3962.t1;Parent=brakerONT.jg3962
#HiC-scaffold_001	AUGUSTUS	start_codon	2844875	2844877	.	+	0	Parent=brakerONT.jg3962.t1
#HiC-scaffold_001	AUGUSTUS	CDS	2844875	2846140	1	+	0	ID=brakerONT.jg3962.t1.cds;Parent=brakerONT.jg3962.t1
#HiC-scaffold_001	AUGUSTUS	stop_codon	2846138	2846140	.	+	0	Parent=brakerONT.jg3962.t1

	tx_by_scaf = defaultdict(list) # key is scaffold, value is list of tx in order
	cds_by_tx = defaultdict(list) # list of CDS bounds by tx
	tx_strand = {}

	# begin parsing
	for line in opentype(brakerfile,'rt'):
		line = line.strip()
		if line and line[0]!="#":
			lsplits = line.split("\t")
			scaffold = lsplits[0]
			feature = lsplits[2]
			feature_start = int(lsplits[3])
			feature_end = int(lsplits[4])
			strand = lsplits[6]
			attributes = lsplits[8]
			attrd = dict([(field.strip().split("=")) for field in attributes.split(";") if field.count("=")])
		#	feature_counter[feature] += 1
			if feature=="gene":
				gene_id = attrd.get("ID")
				gene_buffer[gene_id].append(line)
				gene_order.append(gene_id)
			elif feature=="mRNA":
				gene_id = attrd.get("Parent")
				gene_buffer[gene_id].append(line)
				tx_id = attrd.get("ID")
				tx_strand[tx_id] = strand
				tx_by_scaf[scaffold].append(tx_id)
			elif feature=="CDS":
				tx_id = attrd.get("Parent")
				cds_bounds = (feature_start, feature_end)
				cds_by_tx[tx_id].append(cds_bounds)

				gene_id = tx_id.rsplit(".",1)[0]
				gene_buffer[gene_id].append(line)
			elif feature=="start_codon" or feature=="stop_codon":
				gene_id = attrd.get("Parent").rsplit(".",1)[0]
				gene_buffer[gene_id].append(line)
			else: # should never happen
				pass
	sys.stderr.write("# Collected CDS for {} transcripts\n".format( len(cds_by_tx) ) )

	mRNA_colors = {}
	mRNA_matches = defaultdict(list)

	color_flag_counts = defaultdict(int)

	# green R/F #74c476 #005a32
	# yellow R/F #ffeda0 #fec44f
	# red R/F #fb6a4a #99000d
	sys.stderr.write("# Sorting CDS overlaps to long reads\n")
	for scaffold in sorted(tx_by_scaf.keys()):
		for tx in tx_by_scaf.get(scaffold):
			color_flag = 0
			cds_list = sorted(cds_by_tx.get(tx))
			cds_count = len(cds_list)
			cds_min = min(cb[0] for cb in cds_list)
			cds_max = max(cb[1] for cb in cds_list)
			for lr_bounds, lr_tx in lr_scaf_tx_dict.get(scaffold, {} ).items():
				if lr_bounds[0] > cds_max or lr_bounds[1] < cds_min:
					continue # out of bounds, not relevant for this tx
				lr_exons = lr_tx_exon_dict.get(lr_tx)
				lr_exon_cnt = len(lr_exons)
				i_to_j_overlap = defaultdict(list)
				j_to_i_overlap = defaultdict(list)
				i_to_j_match = defaultdict(list)
				j_to_i_match = defaultdict(list)
				for i,cds in enumerate(cds_list): # [(2851736, 2852444), (2852928, 2853194), (2853477, 2854042)]
					for j,lr_exon in enumerate(sorted(lr_exons)):
						if cds[0] > lr_exon[1]:
							continue
						if cds[1] < lr_exon[0]:
							continue
						if cds[0]==lr_exon[0] or cds[1]==lr_exon[1]: # bound must match exactly
							i_to_j_match[i].append(j)
							j_to_i_match[j].append(i)
						if (lr_exon[0] < cds[0] < lr_exon[1]) or (lr_exon[0] < cds[1] < lr_exon[1]):
							i_to_j_overlap[i].append(j)
							j_to_i_overlap[j].append(i)

				i_to_j_one_one = all_cds_to_one_exon(cds_count , i_to_j_match)
				j_to_i_one_one = all_cds_to_one_exon(lr_exon_cnt, j_to_i_match)
				if i_to_j_one_one and j_to_i_one_one: # meaning matched to any transcript
					color_flag = 1
				if cds_count==1 and all_cds_to_one_exon(cds_count, i_to_j_overlap): # for single exon tx
					color_flag = 1

			color_flag_counts[color_flag] += 1
			if color_flag==0:
				if tx_strand[tx]=="+":
					mRNA_colors[tx] = "DarkRed"
				else:
					mRNA_colors[tx] = "FireBrick"
			elif color_flag==1:
				if tx_strand[tx]=="+":
					mRNA_colors[tx] = "DarkGreen"
				else:
					mRNA_colors[tx] = "DarkSeaGreen"
			else: # should never happen
				pass
	sys.stderr.write("# {} matched, {} unmatched, {} total mRNAs\n".format( color_flag_counts.get(1), color_flag_counts.get(0), len(mRNA_colors) ) )

	# then reprint all genes with color code
	for gene in gene_order:
		for line in gene_buffer.get(gene):
			lsplits = line.split("\t")
			feature = lsplits[2]
			attributes = lsplits[8]
			attrd = dict([(field.strip().split("=")) for field in attributes.split(";") if field.count("=")])
			if feature=="mRNA":
				tx_id = attrd.get("ID")
				new_attributes = "{};color_code={}".format(attributes, mRNA_colors.get(tx_id))
				lsplits[8] = new_attributes
			outline = "{}\n".format( "\t".join(lsplits) )
			sys.stdout.write(outline)

def all_cds_to_one_exon(cds_count, overlaps):
	'''check if all CDS have only 1 hit, return bool'''
	for i in range(cds_count):
		if len(overlaps.get(i,[]))==1:
			continue
		else:
			return False
	else:
		return True

def parse_long_read_tx(longread_gff_file):
	'''parse long read transcripts, return two dicts'''

	tx_by_bounds_by_scaf = defaultdict(dict) # key is scaf, then key is bounds, then value is tx ID
	exons_by_tx = defaultdict(list) # key is tx ID, value is list of exon bounds
	feature_counter = defaultdict(int) # key is feature, value is count

	# HiC-scaffold_001	pinfish	transcript	2423	2615	16	+	.	ID=PIN_pipeline_transcript_60d8e8b8-39ec-4eba-9daa-e740cdb2675d;Name=PIN_pipeline_transcript_60d8e8b8-39ec-4eba-9daa-e740cdb2675d
	# HiC-scaffold_001	pinfish	exon	2423	2615	16	+	.	Parent=PIN_pipeline_transcript_60d8e8b8-39ec-4eba-9daa-e740cdb2675d

	if longread_gff_file.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Reading GFF file {} as gzipped\n".format(longread_gff_file) )
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Reading GFF file {}\n".format(longread_gff_file) )

	# begin parsing
	for line in opentype(longread_gff_file,'rt'):
		line = line.strip()
		if line and line[0]!="#":
			lsplits = line.split("\t")
			scaffold = lsplits[0]
			feature = lsplits[2]
			feature_start = int(lsplits[3])
			feature_end = int(lsplits[4])
			strand = lsplits[6]
			attributes = lsplits[8]
			attrd = dict([(field.strip().split("=")) for field in attributes.split(";") if field.count("=")])

			feature_counter[feature] += 1

			if feature=="transcript": # meaning is top level, should take strand and exons
				gene_id = attrd.get("ID")
				tx_bounds = (feature_start, feature_end)
				tx_by_bounds_by_scaf[scaffold][tx_bounds] = gene_id
			elif feature=="exon":
				gene_id = attrd.get("Parent")
				exon_bounds = (feature_start, feature_end)
				exons_by_tx[gene_id].append(exon_bounds)
	exon_count = sum( len(l) for l in exons_by_tx.values() )
	sys.stderr.write("# File contained:\n")
	for feature, count in sorted(feature_counter.items(), key=lambda x: x[1], reverse=True):
		sys.stderr.write("#{}\t{}\n".format( feature, count ) )
	sys.stderr.write("# Collected {} exons for {} transcripts\n".format( exon_count, len(exons_by_tx) ) )
	return tx_by_bounds_by_scaf, exons_by_tx

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--augustus', help="augustus GFF")
	parser.add_argument('-m','--manual', help="manually annotated sequences, as fasta, can be .gz")
	parser.add_argument('-f','--correspondance-file', help="optional file name to print old names matched with new names")
	parser.add_argument('-l','--gene-label', help="prefix name for ID tag")
	parser.add_argument('-p','--pinfish', help="pinfish GFF")
	parser.add_argument('-s','--stringtie', help="stringtie GFF")
	parser.add_argument('-q','--quiet', action="store_true", help="do not print minor errors")
	args = parser.parse_args(argv)

	main_sc_b_tx_d = defaultdict(dict)
	main_tx_exon_d = defaultdict(list)

	if args.stringtie:
		str_sc_b_tx_d, str_tx_exon_d = parse_long_read_tx(args.stringtie)
		main_sc_b_tx_d.update(str_sc_b_tx_d)
		main_tx_exon_d.update(str_tx_exon_d)

	if args.pinfish:
		pin_sc_b_tx_d, pin_tx_exon_d = parse_long_read_tx(args.pinfish)
		main_sc_b_tx_d.update(pin_sc_b_tx_d)
		main_tx_exon_d.update(pin_tx_exon_d)

	parse_braker(args.augustus, main_sc_b_tx_d, main_tx_exon_d)



if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
