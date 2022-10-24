#!/usr/bin/env python

"""
./fix_final_avas_annotation.py Avas.v1.28_annotations.full_frame_prot.fasta Avas.v1.28_annotations.prot.fasta Avas.v1.28_annotations.gff > Avas.v1.28_annotations.corr.gff
"""

import sys
import time
import re
from collections import defaultdict
from Bio import SeqIO

def get_intervals(intervals, domstart, domlength, doreverse=True):
	'''return a list of intervals with genomic positions for the feature'''
	# example domain arrangement for forward strand
	# intervals from     50,101 127,185 212,300
	# protein domain     71,101 127,185 212,256
	#      in nucleotides  31      59      45
	#      in amino acids  10.3    19.6    15 = 45
	# for domstart at 22 and domlength of 135
	# basestostart is always from transcript N-terminus nucleotide
	# so for forward transcripts, basestostart would be 22, so that 50+22-1=71
	basestostart = int(domstart) # this value always should be 1 or greater
	genomeintervals = [] # will contain a list of tuples
	for interval in sorted(intervals, key=lambda x: x[0], reverse=doreverse):
		intervallength = interval[1]-interval[0]+1 # corrected number of bases
		if basestostart >= intervallength: # ignore intervals before the start of the domain
		#	sys.stderr.write(" ".join([interval, domstart, basestostart, domlength, intervallength]))
			basestostart -= intervallength
		# in example, 101-50+1 = 52, 22 < 52, so else
		else: # bases to start is fewer than length of the interval, meaning domain must start here
			if doreverse: # reverse strand domains
				# if domain continues past an interval, domstart should be equal to interval[1]
				domstart = interval[1] - basestostart + 1 # correct for base numbering at end of interval
				if domstart - interval[0] + 1 >= domlength: # if the remaining part of the domain ends before the start of the interval
					# then define the last boundary and return the interval list
					genomebounds = (domstart-domlength+1, domstart) # subtract remaining length
					genomeintervals.append(genomebounds)
					return genomeintervals
				else:
					genomebounds = (interval[0], domstart)
					genomeintervals.append(genomebounds)
					domlength -= (domstart - interval[0] + 1)
					basestostart = 1 # start at the next interval
			else: # for forward stranded domains
				domstart = interval[0] + basestostart - 1 # correct for base numbering
				if interval[1] - domstart + 1 >= domlength: # if the remaining part of the domain ends before the end of the interval
					# then define the last boundary and return the interval list
					genomebounds = (domstart, domstart+domlength-1) # add remaining length for last interval
					genomeintervals.append(genomebounds)
					return genomeintervals
				else:
					genomebounds = (domstart, interval[1])
					genomeintervals.append(genomebounds)
					domlength -= (interval[1] - domstart + 1)
					basestostart = 1 # next domstart should be interval[0] for next interval
			if domlength < 1: # catch for if all domain length is accounted for
				return genomeintervals
	else:
		sys.stderr.write("WARNING: cannot finish protein at {} for {} in {}\n".format(domstart, domlength, intervals) )
		return genomeintervals




ffprot_file = sys.argv[1]
sys.stderr.write("# reading full frame prots {}\n".format(ffprot_file) )
full_frame_prot_dict = SeqIO.to_dict(SeqIO.parse(ffprot_file,"fasta"))

trimmed_prot_file = sys.argv[2]
sys.stderr.write("# reading trimmed prots {}\n".format(trimmed_prot_file) )
tx_cds_coords = {}
for seqrec in SeqIO.parse(trimmed_prot_file,"fasta"):
	short_prot = str(seqrec.seq)
	short_len = len(short_prot)
	seqid = seqrec.id
	full_frame_prot_frame = int(full_frame_prot_dict[seqrec.id].description.split(" ")[1])
	full_frame_prot = str(full_frame_prot_dict[seqrec.id].seq)
	full_frame_index = full_frame_prot.find(short_prot)
	cds_start = (full_frame_index*3) + 1 + full_frame_prot_frame
	cds_end = cds_start + (short_len)*3 - 1
	try:
		if full_frame_prot[full_frame_index+short_len]=="*":
			cds_end = cds_end+3
	except IndexError:
		pass
	#print( "{}	{}	{}".format(seqid, cds_start, cds_end), file=sys.stdout)
	tx_cds_coords[seqid] = [cds_start, cds_end]

geneintervals = defaultdict(list)
gene_to_strand_dict = {} # key is ID, value is strand as str
gene_to_scaffold_dict = {} # 

gene_order = []

commentlines = 0 # comment lines
linecounter = 0 # all lines that are not comments, even if ignored later
transcounter = 0 # counter for transcript or mRNA
exoncounter = 0 # counter for exon or CDS
ignoredfeatures = 0 # all other features that get ignored

allowed_features = ["gene", "rRNA", "transcript", "exon"]

gtffile = sys.argv[3]

gene_to_tx = defaultdict(list)
lines_by_id = {}
exons_by_tx = defaultdict(list)

programname = ""
sys.stderr.write("# reading GFF {}\n".format(gtffile) )
for line in open(gtffile,'r'):
	line = line.strip()
	if line: # ignore empty lines
		if line[0]=="#": # count comment lines, just in case
			commentlines += 1
		else:
			linecounter += 1
			lsplits = line.split("\t")
			scaffold = lsplits[0]
			programname = lsplits[1]
			feature = lsplits[2]
			strand = lsplits[6]
			attributes = lsplits[8]

			if feature not in allowed_features: # any other features may cause problems later
				ignoredfeatures += 1
				continue

			if attributes.find("ID")>-1: # indicates gff3 format
				geneid = re.search('ID=([\w.|-]+)', attributes).group(1)
			else:
				geneid = None

			if attributes.find("Parent")>-1: # gff3 format but no ID
				toplevel_ID = re.search('Parent=([\w.|-]+)', attributes).group(1)
			else:
				toplevel_ID = None

			if geneid is None and toplevel_ID is None: # means no feature info was found, so error
				raise KeyError("ERROR: cannot extract ID or Parent from line {}\n{}\n".format(linecounter, line) )
			# if either geneid or toplevel_ID are missing, use the other
			if geneid is None and toplevel_ID is not None:
				geneid = toplevel_ID
			if geneid is not None and toplevel_ID is None:
				toplevel_ID = geneid

			if feature=="gene":
				lines_by_id[geneid] = line
				gene_order.append(geneid)

			if feature=="transcript" or feature=="rRNA":
				transcounter += 1
				gene_to_strand_dict[geneid] = strand
				gene_to_scaffold_dict[geneid] = scaffold

				gene_to_tx[toplevel_ID].append(geneid)
				lines_by_id[geneid] = line.replace("transcript","mRNA")

			elif feature=="exon":
				exoncounter += 1
				boundaries = ( int(lsplits[3]), int(lsplits[4]) )
				gene_to_strand_dict[toplevel_ID] = strand
				gene_to_scaffold_dict[toplevel_ID] = scaffold
				geneintervals[toplevel_ID].append(boundaries)
				exons_by_tx[toplevel_ID].append(line)

sys.stderr.write("# Counted {} lines and {} comments  {}\n".format(linecounter, commentlines, time.asctime() ) )
sys.stderr.write("# Counted {} tx and {} exons\n".format( transcounter, exoncounter ) )

querynamedict = defaultdict(int) # counter of unique queries
# count results to filter
not_found_subjects = 0 # counter for subject IDs not found by lookup
total_kept = 0
# count frequency of other problems
missingscaffolds = 0 # count if scaffold cannot be found, suggesting naming problem
intervalproblems = 0 # counter if no intervals are found for some sequence
duplicateintervals = 0 # counter if any queries have duplicate intervals
maxremovals = 0 # counter for hits above max for each query
# count other general stats
intervalcounts = 0
backframecounts = 0
hitDictCounter = defaultdict(int)
linecounter = 0

outputtype = "CDS"

# iterate
sys.stderr.write("# getting CDS coords from exons  {}\n".format( time.asctime() ) )
print( "##gff-version 3", file=sys.stdout )
for gene in gene_order:
	tx_list = gene_to_tx.get(gene)
	print( lines_by_id.get(gene) , file=sys.stdout )
	for tx in tx_list:
		print( lines_by_id.get( tx ) , file=sys.stdout )
		print( "\n".join( exons_by_tx.get(tx) ) , file=sys.stdout )

		hitstart = tx_cds_coords.get(tx,[None,None])[0]
		hitend = tx_cds_coords.get(tx,[None,None])[1]
		if hitstart is None or hitend is None: # has no CDS, do not try to write one
			continue

		hitlength = abs(hitend - hitstart) + 1 # bases 1 to 6 should have length 6
		scaffold = gene_to_scaffold_dict.get(tx, None)
		if scaffold is None:
			missingscaffolds += 1
			if missingscaffolds < 10:
				sys.stderr.write("WARNING: cannot get scaffold for {}\n".format( tx ) )
			elif missingscaffolds == 10:
				sys.stderr.write("WARNING: cannot get scaffold for {}, will not print further warnings\n".format( tx ) )
			continue
		strand = gene_to_strand_dict.get(tx, None)
		genomeintervals = [] # to have empty iterable

		# convert transcript nucleotide to genomic nucleotide, and split at exon bounds
		if strand=='+':
			genomeintervals = get_intervals(geneintervals[tx], hitstart, hitlength, doreverse = False )
		elif strand=='-': # implies '-'
			genomeintervals = get_intervals(geneintervals[tx], hitstart, hitlength, doreverse = True )
		elif strand=='.': # no strand is given by the input GFF
			sys.stderr.write("WARNING: strand is undefined . for {} on {}\n".format(tx, scaffold) )
			continue
		else: # strand is None
			# strand could not be found
			# meaning mismatch between query ID in blast and query ID in the GFF
			sys.stderr.write("WARNING: possible mismatch in ID for {} on {}\n".format(tx, scaffold) )
			continue

		intervalcounts += len(genomeintervals)
		if not len(genomeintervals):
			sys.stderr.write("WARNING: no intervals for {} in {}\n".format(sseqid, tx) )
			intervalproblems += 1
			continue

		# make child features for each interval
		cds_lines = []
		current_phase = 0
		for interval in genomeintervals:
		# thus ID appears as sseqid.tx.number, so avGFP.Renre1234.1, and uses ID in most browsers
			outlist = [scaffold, programname, outputtype, interval[0], interval[1], 1, strand, None, tx]
			cds_lines.append(outlist)


        # phase formula should be (3+p-(l%3))%3
		if strand == "+":
			for fixed_cds_line in sorted(cds_lines, key=lambda x: x[3]):
				fixed_cds_line[7] = current_phase
				#outline = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\tID={8}.cds;Parent={8}\n".format( *fixed_cds_line )
				outline = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\tParent={8}\n".format( *fixed_cds_line )
				sys.stdout.write( outline )
				current_phase = ( 3 + current_phase - (( fixed_cds_line[4]-fixed_cds_line[3] + 1 ) % 3) ) % 3 
		elif strand == "-":
			line_holder = []
			for fixed_cds_line in sorted(cds_lines, key=lambda x: x[4], reverse=True):
				fixed_cds_line[7] = current_phase
				#outline = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\tID={8}.cds;Parent={8}\n".format( *fixed_cds_line )
				outline = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\tParent={8}\n".format( *fixed_cds_line )
				line_holder.append( outline )
				current_phase = ( 3 + current_phase - (( fixed_cds_line[4]-fixed_cds_line[3] + 1 ) % 3) ) % 3 
			for outline in line_holder[::-1]:
				sys.stdout.write( outline )
		else:
			sys.stderr.write("WARNING: strand is undefined . for {} on {}\n".format(tx, scaffold) )


sys.stderr.write("# found {} CDS intervals\n".format( intervalcounts ) )

if missingscaffolds:
	sys.stderr.write("# WARNING: could not find scaffold for {} hits  {}\n".format(missingscaffolds, time.asctime() ) )
if intervalproblems:
	sys.stderr.write("# WARNING: {} matches have hits extending beyond gene bounds  {}\n".format(intervalproblems, time.asctime() ) )




