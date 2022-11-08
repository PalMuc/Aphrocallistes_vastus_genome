#!/usr/bin/env python
#
# gtfstats.py v1.0 created 2015-07-27
# original version made for Francis and Woerheide 2017 GBE
# 
# current version available at
# https://bitbucket.org/wrf/sequences/src/master/

'''
gtfstats.py v1.8  last modified 2022-03-21
    get summary information about annotation from GTF or GFF

gtfstats.py -i genes.gtf

    to additionally calculate gap positions relative to assembly size, use -s

gtfstats.py -i genes.gtf -s scaffolds.fasta

    if gtf is not in format of gene/exon or transcript/exon
    use -g option to calculate genes based on exon groups
    note that this may not work, depending on GFF format

    if exons are not specified, such as for intron/CDS format
    use -c option to calculate exons from CDS
    and -u to also include five and three prime UTR

    use -m to include match or cDNA_match types as exons, such as
    from transcripts mapped to a genome

    if both gene and mRNA are specified, use -G to ignore gene types
    this is typical for Ensembl GFF files
    otherwise makes this error:
  geneid = re.search('[\w_]+ "([\w.|-]+)";', attributes).group(1)
  AttributeError: 'NoneType' object has no attribute 'group'

    for mapped Trinity genes, use -T to split c1234_g1_i1.path1 at c1234_g1

    for AUGUSTUS, use -c and -G since exons are not specified

    most UnboundLocalErrors are solved by one of the above options
'''

#
import sys
import argparse
import time
import re
import os
import gzip
from collections import defaultdict,Counter
from Bio import SeqIO
#

def combine_intervals(rangelist):
	'''convert list of tuples to non redundant invervals'''
	# sl = [(1,10), (1,6), (15,20), (1,10), (19,29), (30,35), (6,13), (40,48), (15,21), (42,50)]
	nrintervallist = []
	srtrangelist = sorted(rangelist) # sort list now, to only do this once
	interval = srtrangelist[0] # need to start with first time, which will be the same for the first bounds
	for bounds in srtrangelist:
		# since it is sorted bounds[0] should always be >= interval[0]
		if bounds[0] > interval[1]+1: # if the next interval starts past the end of the first + 1
			nrintervallist.append(interval) # add to the nr list, and continue with the next
			interval = bounds
		else:
			if bounds[1] > interval[1]: # bounds[1] <= interval[1] means do not extend
				interval = (interval[0], bounds[1]) # otherwise extend the interval
	else: # append last interval
		nrintervallist.append(interval)
	# should return [(1, 13), (15, 35), (40, 50)]
	return nrintervallist


def sum_interval_span(rangelist):
	"""from list of tuples, return sum of all intervals"""
	interval_sum = 0
	for bounds in rangelist:
		interval_sum = bounds[1] - bounds[0] + 1
	return interval_sum


def get_occupied_length(rangedict): # added in v1.3 for speed and memory improvements
	'''from a dictionary of list, where list items are tuples of intervals, return the sum of non redundant intervals'''
	# rangedict must be dictionary of lists
	# key is scaffold
	# where values are list items are boundaries of features
	sys.stderr.write("# Determining occupied bases  {}\n".format( time.asctime() ) )
	occlen = 0
	occintervals = 0
	for rl in rangedict.values(): # rl in rangedict.values() should be as [(1,5), (25,29)]
		for nrintvl in combine_intervals(rl):
			occintervals += 1
			occlen += nrintvl[1]+1-nrintvl[0] # calculate adjusted interval length
	sys.stderr.write("Counted {} non-redundant groups (when merging)\n".format(occintervals) )
	return occlen


def get_exon_stats(exonboundaries, genekey, print_exon_counts, print_exons=False):
	'''given a dict of lists of tuple intervals, count the non redundant sum of intervals'''
	exoncountbygroup = [] # list of integers, number of exons per gene/transcript
	for geneid, exonlist in exonboundaries.items():
		exonset = set(exonlist) # removes redundant exons when counting by gene
		exoncountbygroup.append(len(exonset))
		if print_exons: # display as g3773.t1 8581960 8583538 1579
			for exon in sorted(exonset):
				sys.stdout.write( "{} {} {} {}\n".format(geneid, exon[0], exon[1], exon[1]-exon[0]+1) )
	totalexons = sum(exoncountbygroup)
	avgexnum = totalexons * 1.0 / len(exoncountbygroup)
	sys.stderr.write("Counted {} non-redundant exons (unmerged at boundaries)\n".format(totalexons) )
	sys.stderr.write("Average {:.6f} exons per {}, max {}\n".format(avgexnum, genekey, max(exoncountbygroup) ) )
	if print_exon_counts: # print the counter of exons per gene
		sys.stderr.write("{}".format(Counter(exoncountbygroup)) )
	return totalexons



def get_intron_stats(nametoscaffold, exonboundaries, do_print_intervals, get_nested, interval_desc):
	'''from exon intervals by gene, return bulk parameters about introns or intergenic, and/or nested genes'''
	introncount = 0
	intronlengths = []

	is_introns = ( interval_desc=="introns" )

	introns_by_gene = defaultdict(list)
	introns_by_scaffold = defaultdict(list) # key is scaffold, value is list of intron intervals on that scaffold
	genebounds_by_name = {} # key is gene name, value is bounds, needed for exon counting later
	exonic_bp_by_gene = {} # key is gene name, value is integer of summed non-redundant exons

	sys.stderr.write("# Measuring {} stats  {}\n".format( interval_desc, time.asctime() ) )
	for geneid, exonlist in exonboundaries.items():
		exonset = set(exonlist) # removes redundant exons when counting by gene
		newbounds = (min(ex[0] for ex in exonlist), max(ex[1] for ex in exonlist) )
		genebounds_by_name[geneid] = newbounds
		nr_exonset = combine_intervals( exonlist )
		exonic_length = sum_interval_span(nr_exonset)
		exonic_bp_by_gene[geneid] = exonic_length
		if len(exonset) > 1:
			scaffold = nametoscaffold.get(geneid, None)
			introncount += len(exonset)-1 # number of introns should be exons minus 1 per transcript
			sortedexonlist = sorted(exonset) # should by default sort first number then second, no need for key=lambda x: x[0]
			# iterate through pairs of intervals
			for x,y in zip(sortedexonlist[:-1],sortedexonlist[1:]):
				# from exons (1,5) and (25,29), intron should be (6,24)
				# so length is 25-5-1, or 24-6+1
				intron_start = x[1]+1
				intron_end = y[0]-1
				intronlen = intron_end - intron_start + 1 # this works regardless of transcript direction
				# exons are sorted, so checks for negative intron length
				# will only be less than 1 if gene has two exons with different 3-prime splice sites
				# like (10,15) and (10,20), tries to calculate 10-15
				# hence considers bases (16,20) as exon, not intron
				if intronlen > 0:
					intronlengths.append(intronlen)
					introns_by_gene[geneid].append( (intron_start, intron_end) )
					introns_by_scaffold[scaffold].append( (intron_start, intron_end) )
					if do_print_intervals: # display as transcript:FBtr0078764 (5493105, 5493489) (5498777, 5498870) 5287
						sys.stdout.write( "{} {} {} {}\n".format(geneid, x, y, intronlen) )
	# sum list of values of lengths
	intronsum = sum(intronlengths)

	# for cases where there are no introns, and divide by zero fails
	if not intronsum: # for example, on a bacterial GFF
		return 0,0,0

	avgintron = intronsum * 1.0 / len(intronlengths)
	sys.stderr.write("Counted {} {}, added {}\n".format(introncount, interval_desc, len(intronlengths) ) )

	# calculation assumes introns were given
	if is_introns and get_nested:
		# this step differs from the analysis by
		# https://github.com/conchoecia/chep/blob/master/scripts/gff_to_intron_bed.py
		# which marginally trims the exons to increase overlap, and
		# removes a fraction of the longest introns, 
		# which may be due to trans-spliced leaders remapping to the genome
		sys.stderr.write("# Counting nested genes  {}\n".format( time.asctime() ) )
		nested_gene_count = 0
		nested_gene_bases = 0
		nested_exon_bases = 0

		# for each gene, determine if it fits completely within one or more introns
		# to prevent any transcripts from being double counted
		# THIS STEP IS SLOW #TODO
		for geneid, geneinterval in genebounds_by_name.items():
			scaffold = nametoscaffold.get(geneid, None)
			is_nested = 0
			# check against all introns on that scaffold
			for introninterval in introns_by_scaffold.get(scaffold,[]):
				if geneinterval[0] <= introninterval[0]:
					continue
				if geneinterval[1] >= introninterval[1]:
					continue
				is_nested += 1
			# if was flagged as nested one or more times
			if is_nested > 0: # can be more than one #TODO
				#print("{}\t{}\t{}\t{}\t{}".format(scaffold, geneid, is_nested, geneinterval[0], geneinterval[1]), file=sys.stderr)
				gene_length = geneinterval[1] - geneinterval[0] + 1
				nested_gene_count += 1
				nested_gene_bases += gene_length
				nested_exon_bases += exonic_bp_by_gene.get(geneid)
		sys.stderr.write("Counted {} nested genes, of {} total bases, {} exonic bases  {}\n".format(nested_gene_count, nested_gene_bases, nested_exon_bases, time.asctime() ) )

	return intronsum, avgintron, len(intronlengths)


def exons_to_transcript_bounds(nametoscaffold, exonboundaries):
	'''from exon intervals by gene and scaffold by exons, return a dict where keys are scaffolds and values are lists of gene intervals'''
	sys.stderr.write("# Determining genes from exons  {}\n".format( time.asctime() ) )
	genelengths = [] # list of lengths
	genesbyscaffold = defaultdict(list) # keys are scaffold names, values are lists of tuples of boundaries
	for geneid, exonlist in exonboundaries.items():
		# from list of exon bounds, like (150,250),(500,700) should return (150,700)
		newbounds = (min(ex[0] for ex in exonlist), max(ex[1] for ex in exonlist) )
		genelen = (newbounds[1] - newbounds[0] + 1 )
		genelengths.append(genelen)
		scaffold = nametoscaffold[geneid]
		genesbyscaffold[scaffold].append(newbounds)
	return genelengths, genesbyscaffold


def report_scaffold_stats(genesbyscaffold):
	'''across all scaffolds, get counts for average and median number of genes'''
	num_GpS = [] # number of genes per scaffold
	max_GpS = 0
	max_GpS_scaf = ""
	for scaffold, intervals in genesbyscaffold.items():
		GpS = len(intervals)
		num_GpS.append(GpS)
		if GpS > max_GpS:
			max_GpS = GpS
			max_GpS_scaf = scaffold
	total_genes = sum(num_GpS)
	total_scaffolds_w_genes = len(genesbyscaffold)
	average_genes = total_genes * 1.0 / len(num_GpS)
	median_genes = sorted(num_GpS)[len(num_GpS)//2]
	sys.stderr.write("{} total genes on {} scaffolds, max was {} on {}  \n".format(total_genes, total_scaffolds_w_genes, max_GpS, max_GpS_scaf ) )
	sys.stderr.write("mean {:.2f} genes per scaffold, median {} genes per scaffold  \n".format(average_genes, median_genes) )
	return total_scaffolds_w_genes


def attributes_to_dict(attributes):
	'''convert GFF attribute string into dictionary of key-value pairs'''
	attrd = {}
	if attributes.find("ID=")>-1 or attributes.find("Parent=")>-1: # indicates GFF3 format
		# if one of the terms does not have = sign, perhaps Note, then ignore
		attrd = dict([(field.strip().split("=",1)) for field in attributes.split(";") if field.count("=")])
	else: # assume GTF format
		try:
			attrd = dict([(field.strip().split(" ",1)) for field in attributes.split(";") if field])
		except ValueError: # catch for Ensembl genomes, which use = but not ID
			attrlist = [field for field in attributes.split(";") if field]
			for attr in attrlist:
				try:
					if attr.count("=")>0:
						attrd.update(dict([attr.strip().split("=")]))
					elif attr.count(" ")>0:
						attrd.update(dict([attr.strip().split(" ")]))
					else: # apparently the field is not delimited
						attrd["NULL"] = attr
				except ValueError: # for in line comments like some Broad Institute gtfs
					# '# At least one base has a quality score < 10'
					sys.stderr.write("WARNING: UNKNOWN ATTRIBUTE: {}\n".format(attr) )
	return attrd


def find_gaps(scaffoldfile, exclusiondict, repeatletter, genesbyscaffold, get_id_from_description, verbose, stranded, above=2, below=1000000000):
	'''read in fasta scaffolds and determine if gaps are within boundaries of genes'''
	repeatregex = re.compile("([{0}{1}])+".format(repeatletter, repeatletter.lower() ) )

	sys.stderr.write("# Determining gap positions from scaffolds {}  {}\n".format( scaffoldfile, time.asctime() ) )

	if get_id_from_description:
		sys.stderr.write("# Extracting scaffold IDs from fasta description, instead of default ID\n")

	scaflength = 0 # total length of all scaffolds combined
	contigcounter = 0
	repcounter = 0
	geneoverlaps = 0
	introngaps = 0
	intergenicgaps = 0
	length_list = [] # list of scaffold lengths, for n50 calculation

	if scaffoldfile.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Reading gaps, opening {} as gzipped  {}\n".format(scaffoldfile, time.asctime() ) )
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Reading gaps, opening {}  {}\n".format(scaffoldfile, time.asctime() ) )
	# iterate over sequences
	for seqrec in SeqIO.parse(opentype(scaffoldfile,'rt'), "fasta"):
		contigcounter += 1
		seqlength = len(seqrec.seq)
		scaflength += seqlength
		length_list.append(seqlength)
		if get_id_from_description: # assembly is from GenBank, original scaffold IDs are in description
			# ID should be in format of:
			#>QEPC01000001.1 Sphenodon punctatus isolate mauimua-1 ScrUdWx_8519, whole genome shotgun sequence
			# original scaffold ID of the submission should be ScrUdWx_8519
			contig = seqrec.description.split(",",1)[0].split(" ")[-1]
		else: # use whatever ID is parsed normally
			contig = seqrec.id

		# check if on the exclusion list
		if exclusiondict and contig in exclusiondict:
			continue

		# stranded only operations
		if stranded: # if stranded, check forward strand, and then reverse
			contig += "+" # first contig is plus strand
			scaflength += len(seqrec.seq) # length must be double counted

		# resume normal operation
		genelistbyscaffold = genesbyscaffold[contig]
		for rep in repeatregex.finditer(str(seqrec.seq)): # iterate through normal repeats
			replen = rep.end() - rep.start()
			if below > replen >= above:
				repcounter += 1
				repspan = [x+1 for x in rep.span()] # correct python index to GFF
				for genebound in genelistbyscaffold:
					if repspan[0] > genebound[0] and repspan[1] < genebound[1]:
						introngaps += replen
						break # if a gene is found, stop searching
					elif (repspan[0] < genebound[0] and repspan[1] > genebound[0]) or (repspan[0] < genebound[1] and repspan[1] > genebound[1]):
						if verbose:
							sys.stderr.write("WARNING: gap at {} spans gene {} on {}\n".format(repspan, genebound, contig) )
						geneoverlaps += 1
						intergenicgaps += replen # these are always at boundaries of genes by 1 or 2 bases
						break
				else: # is therefore otherwise in an intergenic region
					intergenicgaps += replen

				# block is repeat of above block, but for minus strand in stranded mode
				if stranded:
					for genebound in genesbyscaffold[contig+"-"]:
						if repspan[0] > genebound[0] and repspan[1] < genebound[1]:
							introngaps += replen
							break
						elif (repspan[0] < genebound[0] and repspan[1] > genebound[0]) or (repspan[0] < genebound[1] and repspan[1] > genebound[1]):
							if verbose:
								sys.stderr.write("WARNING: gap at {} spans gene {} on {}\n".format(repspan, genebound, contig) )
							geneoverlaps += 1
							intergenicgaps += replen
							break
					else:
						intergenicgaps += replen
				# end of stranded search

	sys.stderr.write("Counted {} total bases for {} contigs  {}\n".format( scaflength, contigcounter, time.asctime() ) )
	sys.stderr.write("Found {} bases of gaps within genes (introns)\n".format(introngaps) )
	sys.stderr.write("Found {} bases of gaps between genes (intergenic)\n".format(intergenicgaps) )
	if geneoverlaps:
		sys.stderr.write("Found {} gaps overlapping with gene boundaries\n".format(geneoverlaps) )

	if get_id_from_description:
		sys.stderr.write("# Scaffold IDs extracted as {} from {}\n".format(contig, seqrec.description) )

	# calculate scaffold N50
	halfsum = scaflength/2
	ncountsum = 0
	scaffold_n50 = 0
	for i_len in sorted(length_list):
		# changed to output float (%.2f, not %d) for v2.3
		ncountsum += i_len
		if ncountsum >= halfsum:
			scaffold_n50 = i_len
			break

	# send back sclen, igaps, ggaps
	return contigcounter, scaffold_n50, scaflength, introngaps, intergenicgaps


def find_repeats( repeatmaskergff, genesbyscaffold, exclusiondict ):
	'''read GFF of repeat masker output, and return two dicts, where keys are repeat classes, and values are total bases of each class'''
	intron_repeats = defaultdict(int) # count of number of repeats
	intron_bases = defaultdict(int) # total bases of repeats found in introns
	intergenic_repeats = defaultdict(int)
	intergenic_bases = defaultdict(int)
	repcounter = 0
	geneoverlaps = 0
	if repeatmaskergff.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Reading repeats from {} as gzipped  {}\n".format(repeatmaskergff, time.asctime() ) )
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Reading repeats from {}  {}\n".format(repeatmaskergff, time.asctime() ) )

	# example from Carcharodon genome
	#ScRDnwc_29106_354205	repeatmasker	match	2067	2097	236	+	.	ID=ScRDnwc_29106_354205:hit:48490:1.3.0.0;Name=species:SacSINE1|genus:SINE%2FtRNA-Deu-L2;Target=species:SacSINE1|genus:SINE%2FtRNA-Deu-L2 108 138 +
	#ScRDnwc_29106_354205	repeatmasker	match_part	2067	2097	236	+	.ID=ScRDnwc_29106_354205:hsp:70186:1.3.0.0;Parent=ScRDnwc_29106_354205:hit:48490:1.3.0.0;Target=species:SacSINE1|genus:SINE%252FtRNA-Deu-L2 108 138 +

	# read only match features with program repeatmasker
	for line in opentype(repeatmaskergff,'r'):
		line = line.strip()
		if line and line[0]!="#": # ignore empty lines and comments
			lsplits = line.split("\t")
			if len(lsplits) < 9: # random other weird lines in the gff
				continue
			scaffold = lsplits[0]
			if exclusiondict and exclusiondict.get(scaffold, False):
				continue # skip anything that hits to excludable scaffolds
			program = lsplits[1]
			feature = lsplits[2]
			if program=="repeatmasker" and feature=="match":
				repcounter += 1
				repstart = int(lsplits[3])
				repend = int(lsplits[4])
				replen = repend - repstart + 1
				attributes = lsplits[8]
				attrd = dict([(field.strip().split("=")) for field in attributes.split(";") if field.count("=")])
				repgenus = attrd["Name"].split(":")[-1].replace("%2F","-")
				genelistbyscaffold = genesbyscaffold[scaffold]
				for genebound in genelistbyscaffold:
					if repstart > genebound[0] and repend < genebound[1]:
						intron_repeats[repgenus] += 1
						intron_bases[repgenus] += replen
						break # if a gene is found, stop searching
					elif (repstart < genebound[0] and repend > genebound[0]) or (repstart < genebound[1] and repend > genebound[1]):
						geneoverlaps += 1
						intergenic_repeats[repgenus] += 1
						intergenic_bases[repgenus] += replen # probably at boundaries of genes by 1 or 2 bases
						break
				else: # not clearly in a gene, so is therefore in an intergenic region
					intergenic_repeats[repgenus] += 1
					intergenic_bases[repgenus] += replen
				if not repcounter % 2000:
					sys.stderr.write(".")
				if not repcounter % 100000:
					sys.stderr.write(repcounter + time.asctime() + os.linesep)
	else: # print final count
		sys.stderr.write(repcounter + time.asctime() + os.linesep)
	sys.stderr.write("Counted {} total repeats for {} bases\n".format(repcounter, sum(intron_bases.values()) + sum(intergenic_bases.values()) ) )
	if geneoverlaps:
		sys.stderr.write("Found {} gaps overlapping with gene boundaries\n".format(geneoverlaps) )
	# make dict of unique repeat types
	repeat_types = {}
	for reptype in ( intron_repeats.keys() + intergenic_repeats.keys() ):
		repeat_types[reptype] = True
	sys.stderr.write("Found {} repeat types\n".format(len(repeat_types)) )
	for reptype in repeat_types.keys():
		sys.stdout.write("{}\t{}\t{}\t{}\t{}\n".format(reptype, intron_repeats.get(reptype,0), intron_bases.get(reptype,0), intergenic_repeats.get(reptype,0), intergenic_bases.get(reptype,0) ) )
	# no return


def make_exclude_dict(excludefile):
	'''read file of list of contigs, and return a dict where keys are contig names to exclude'''
	sys.stderr.write("# Reading exclusion list {}  {}\n".format(excludefile, time.asctime() ) )
	exclusion_dict = {}
	for term in open(excludefile,'r'):
		term = term.strip()
		if term[0] == ">":
			term = term[1:]
		exclusion_dict[term] = True
	sys.stderr.write("# Found {} contigs to exclude  {}\n".format(len(exclusion_dict), time.asctime() ) )
	return exclusion_dict


def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-i','--input', help="input gtf/gff3 file", required=True)
	parser.add_argument('-E','--exclude', help="file of list of bad contigs")
	parser.add_argument('-c','--CDS', action="store_true", help="exons are missing, use CDS instead")
	parser.add_argument('-d','--stranded', action="store_true", help="calculate strands separately")
	parser.add_argument('-g','--no-genes', action="store_true", help="genes are not defined, get gene ID for each exon")
	parser.add_argument('-G','--ignore-gene', action="store_true", help="skip lines where type is gene, take gene information from transcript or mRNA types")
	parser.add_argument('-m','--match', action="store_true", help="treat match and cDNA_match as exons")
	parser.add_argument('-M','--mid-level-attr', default="ID", help="attribute to identify mid-level features, mRNA/transcripts [ID]")
	parser.add_argument('-N','--nested-genes', action="store_true", help="report number and amount of nested genes")
	parser.add_argument('-O','--ignore-introns', action="store_true", help="disregard intron features, and count from exons and genes")
	parser.add_argument('-p','--keep-pseudo', action="store_true", help="do not ignore pseudogenes")
	parser.add_argument('-r','--repeat', default="N", help="gaps counted as N or n")
	parser.add_argument('-s','--scaffolds', help="scaffolds as fasta for gap calculations, can be .gz")
	parser.add_argument('-t','--by-transcript', metavar="TYPE", nargs="?", const="transcript", default="gene", help="count by transcript id instead of gene id")
	parser.add_argument('-T','--trinity', action="store_true", help="gene IDs are from Trinity")
	parser.add_argument('-u','--UTR', action="store_true", help="treat UTRs as exons if exons are missing")
	parser.add_argument('-v','--verbose', action="store_true", help="verbose output")
	parser.add_argument('-w','--w', action="store_true", help="give summarized output")
	parser.add_argument('--exon-counter', action="store_true", help="print histogram data of exons per gene")
	parser.add_argument('--print-exons', action="store_true", help="print exon information to stdout")
	parser.add_argument('--print-introns', action="store_true", help="print intron information to stdout")
	parser.add_argument('--print-intergenic', action="store_true", help="print intergenic distance information to stdout")
	parser.add_argument('--use-original-scaffolds', action="store_true", help="for a GenBank assembly fasta file, extract the scaffold IDs from the fasta description instead of using the normal ID")
	parser.add_argument('--repeat-masker', help="GFF file of repeat masker features, for additional calculations")
	parser.add_argument('--debug', action="store_true", help="debug mode, print various tables")
	args = parser.parse_args(argv)

	commentlines = 0

	featurecounts = defaultdict(int)
	exonsum = 0
	intronsum = 0
	#genesum = 0 # no longer needed with list of lengths
	matchsum = 0
	genelengths = [] # list of lengths, to calculate mean min and max
	mrnalengths = defaultdict(int)

	exonbyscaffold = defaultdict(list) # make list of exon intervals, for the count of exonic sequence
	genesbyscaffold = defaultdict(list) # make list of gene intervals, for the count of genic sequence
	exonboundaries = defaultdict(list) # make list of tuples of exons by transcript, to later calculate introns

	nscaf = 0 # will report number of scaffolds with genes, otherwise will report all scaffolds if fasta is given
	scafn50 = 0 # dummy value, to be reassigned if fasta scaffolds are given

	# in order to get transcript boundaries with -g, store names to scaffolds
	nametoscaffold = {}

	idrekey = args.by_transcript # either "gene" or "transcript"

	exclusiondict = make_exclude_dict(args.exclude) if args.exclude else None

	pseudodict = {} # for pseudogenes, key is geneid, value is True
	geneid = "" # empty string so can be called if no transcripts are given but exons are needed

	transcript_features = ["transcript", "mRNA", "rRNA", "tRNA", "ncRNA", "snoRNA", "lincRNA", "lnc_RNA", "miRNA"]
	utr_features = ["five_prime_UTR", "three_prime_UTR", "UTR", "UTR_3", "UTR_5"] # UTR_3 and UTR_5 types are non-standard

	# indicate some flags
	if args.no_genes:
		sys.stderr.write("# Genes are undefined, searching IDs for each feature\n")
	if args.ignore_gene:
		sys.stderr.write("# Ignoring gene features -G \n")
	if args.ignore_introns:
		sys.stderr.write("# Ignoring intron features -O \n")
	if args.keep_pseudo:
		sys.stderr.write("# Keeping pseudogene features for analysis\n")
	if args.stranded:
		sys.stderr.write("# Counting each strand separately\n")

	### MAIN LOOP ###
	if args.input.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Parsing {} as gzipped  {}\n".format(args.input, time.asctime() ) )
	else: # otherwise assume normal open for GTF format
		opentype = open
		sys.stderr.write("# Parsing {}  {}\n".format(args.input, time.asctime() ) )
	for line in opentype(args.input,'rt'):
		line = line.strip()
		if line: # ignore empty lines
			if line[0]=="#": # count comment lines, just in case
				commentlines += 1
			else:
				lsplits = line.split("\t")
				if len(lsplits) < 9: # random other weird lines in the gff
					continue
				scaffold = lsplits[0]
				if args.exclude and exclusionDict.get(scaffold, False):
					continue # skip anything that hits to excludable scaffolds
				if args.stranded:
					scaffold = "{}{}".format(scaffold, lsplits[6])
				feature = lsplits[2]
				featurecounts[feature] += 1
				if args.verbose: # this only removes comment lines
					wayout.write( line + os.linesep )
				attributes = lsplits[8]

				### GENE TYPE ###
				if feature=="gene" and not args.ignore_gene:
					# canonically this ignores miRNA_gene pseudogene rRNA_gene snRNA_gene snoRNA_gene tRNA_gene ncRNA_gene
					attrd = attributes_to_dict(attributes)
					try:
						#geneid = re.search('ID=([\w.|-]+);', attributes).group(1)
						geneid = attrd["ID"]
					except KeyError: # meaning no ID, so not GFF3 somehow
						geneid = attrd.get("{}_id".format(idrekey),None)
						if geneid==None and feature=="gene":
							raise KeyError("ERROR: {}_id not found, gene features appear unformatted, try to rerun with -G".format(idrekey) )
					if "pseudo" in attrd and attrd["pseudo"] == "true":
						pseudodict[geneid] = True
						if not args.keep_pseudo:
							continue
					if "partial" in attrd and attrd["partial"] == "true":
						pseudodict[geneid] = True
						if not args.keep_pseudo:
							continue
					genelen = int(lsplits[4]) - int(lsplits[3]) + 1
					genelengths.append(genelen)
					boundaries = (int(lsplits[3]), int(lsplits[4]) ) # tuple of ints
					genesbyscaffold[scaffold].append(boundaries)
					gb = boundaries
					gl = geneid

				### GFF DEFINED PSEUDOGENES ###
				elif feature=="pseudogene" or feature=="pseudo_gene":
					attrd = attributes_to_dict(attributes)
					try:
						geneid = attrd["ID"]
					except KeyError: # meaning no ID, so not GFF3 somehow
						geneid = attrd["{}_id".format(idrekey)]
					pseudodict[geneid] = True

				### TRANSCRIPT TYPES ###
				elif feature in transcript_features: # if feature is in the list of allowed transcript types
				# this ignores RNA miRNA snRNA snoRNA pseudogenic_tRNA
					attrd = attributes_to_dict(attributes)
					try:
						#geneid = re.search('ID=([\w.|:-]+);', attributes).group(1) # : added to allow for matching in Octopus gene models
						geneid = attrd[args.mid_level_attr] # default is "ID", could be "Parent
						# this would identify redundant exons between transcripts
					except KeyError: # meaning no ID, so not GFF3 somehow
						if "{}_id".format(idrekey) in attrd:
						# geneid = re.search('{}_id "([\w.|+-]+)";'.format(idrekey), attributes).group(1) # + added to allow for matching in Mnemiopsis cufflinks transcripts, maybe not needed
							geneid = attrd["{}_id".format(idrekey)].replace('"','') # remove " from stringtie
						elif "name" in attrd: # for most other GTF cases, such as name, Name, mrna
							geneid = attrd["name"]
						else: # for some AUGUSTUS versions, just take the straight attributes
							geneid = attributes
					if args.trinity: # for Trinity/GMAP mapping
						geneid = geneid.rsplit("_",1)[0] # split c1234_g1_i1.path1 at c1234_g1
					try:
						parent = attrd["Parent"]
						if pseudodict.get(parent,False) and not args.keep_pseudo:
							pseudodict[geneid] = True
							continue
					except KeyError:
						pass # apparently do nothing
					translen = int(lsplits[4]) - int(lsplits[3]) + 1
					genelengths.append(translen)
					boundaries = (int(lsplits[3]), int(lsplits[4]) ) # tuple of ints
					genesbyscaffold[scaffold].append(boundaries)

				### EXON TYPES AND SUBSTITUTES ###
				elif feature=="exon" or (args.CDS and feature=="CDS") or (args.UTR and feature in utr_features) or (args.match and feature=="cDNA_match"):
					# might need (args.UTR and feature=="5'-UTR") or (args.UTR and feature=="3'-UTR") for Branchiostoma
					attrd = attributes_to_dict(attributes)
					if pseudodict: # if any items are pseudogenes, start checking exons
						try:
							parent = attrd["Parent"]
							if pseudodict.get(parent,False) and not args.keep_pseudo:
								continue # skip exon if the parent is partial or pseudogene
						except KeyError:
							pass # apparently do nothing
					exonlength = int(lsplits[4]) - int(lsplits[3]) + 1
					exonsum += exonlength
					if args.no_genes:
						geneid = attrd.get("ID", None)
						if not geneid:
							try: # for GMAP output format
								geneid = re.search('ID=([\w.|_-]+);', attributes).group(1)
							except AttributeError: # in case re fails and group does not exist
								if "gene_id" in attrd:
									geneid = attrd["gene_id"].replace('"','') # remove " from stringtie
								elif "name" in attrd: # for JGI Emihu1 reduced gene format
									geneid = attrd["name"]
								else:
								# this should work for basically all other types, gene_id and name
									geneid = re.search('[\w_]+ "([\w._|-]+)";', attributes).group(1)
							#geneid = re.search('name "([\w.|-]+)";', attributes).group(1)
					elif not geneid: # geneid is still needed for some calculations
						# take from Parent if somehow a geneid does not carry from a transcript feature
						geneid = attrd.get("Parent", None)
						if geneid is None:
							raise KeyError("WARNING: cannot get gene ID for exon types from transcript type, try with -g")

					nametoscaffold[geneid] = scaffold

					mrnalengths[geneid] += exonlength
					exonbounds = (int(lsplits[3]), int(lsplits[4]))
					exonbyscaffold[scaffold].append(exonbounds) # for occupied exon space for scaffolds
					# exonboundaries cannot be used to get occupied space, since it ignores overlapping genes
					exonboundaries[geneid].append(exonbounds) # for calculating intron boundaries, by gene
				elif feature=="intron":
					intronlength = int(lsplits[4]) - int(lsplits[3]) + 1
					intronsum += intronlength
	sys.stderr.write("Counted {} comment lines  ".format(commentlines) + time.asctime() + os.linesep)
	if pseudodict:
		sys.stderr.write("Flagged {} genes or transcripts as pseudogenes\n".format(len(pseudodict)) )

	### OUTPUT ###
	for k in sorted(featurecounts.keys()):
		sys.stderr.write("{}\t{}\n".format(k, featurecounts[k]) )
	if exonsum:
		if args.CDS: # for gffs with no defined exons, using CDS instead
			avgexon = exonsum * 1.0 / (featurecounts["CDS"]+featurecounts["exon"])
		elif args.match: # for gffs from GMAP where exons are not defined
			avgexon = exonsum * 1.0 / featurecounts["cDNA_match"]
		else: # for normal exon gtfs
			avgexon = exonsum * 1.0 / featurecounts["exon"]

		sys.stderr.write("Total exon length {} , average {:.2f}\n".format(exonsum, avgexon) )
		exonOccupiedLength = get_occupied_length(exonbyscaffold)
		sys.stderr.write("Total non-redundant exon length {}  {}\n".format(exonOccupiedLength, time.asctime() ) )

		nrexoncount = get_exon_stats(exonboundaries, idrekey, args.exon_counter, args.print_exons)
		sys.stderr.write("Average non-redundant exon length {:.2f}  {}\n".format(exonOccupiedLength*1.0/nrexoncount , time.asctime() ) )
	else: # suggest alternative if no exon features are found
		exonOccupiedLength = 0
		if featurecounts["CDS"]:
			sys.stderr.write("# WARNING: exon sum is 0, CDS features detected, re-run with option -c\n")
		if featurecounts["cDNA_match"]:
			sys.stderr.write("# WARNING: exon sum is 0, match features detected, re-run with option -m\n")

	if args.no_genes: # if genes were not specified, get gene info from exons
		genelengths, genesbyscaffold = exons_to_transcript_bounds(nametoscaffold, exonboundaries)

	if intronsum and not args.ignore_introns: # if introns were defined in the gtf, use them
		avgintron = intronsum * 1.0 / featurecounts["intron"]
		introncount = featurecounts["intron"]
		sys.stderr.write("Total intron length from features {} , average {:.2f}\n".format(intronsum, avgintron) )
	elif exonboundaries: # otherwise introns were not defined in the gtf, so calculate from exons
		intronsum, avgintron, introncount = get_intron_stats( nametoscaffold, exonboundaries, args.print_introns, args.nested_genes, "introns" )
		sys.stderr.write("Total intron length inferred from exons {} , average {:.2f}\n".format(intronsum, avgintron) )

	if matchsum: # count of cDNA_match, may not be relevant for properly annotated genomes
		avgexon = matchsum * 1.0 / featurecounts["cDNA_match"]
		sys.stderr.write("Total match length {} , average {:.2f}\n".format(matchsum, avgexon) )

	if mrnalengths:
		mrnacount = len(mrnalengths)
		sys.stderr.write("Counted {} {}s (across all transcript types)\n".format(mrnacount, idrekey) )
		mrnasum = sum(mrnalengths.values())
		longest_mrna = max(mrnalengths.values())
		avgmrna = mrnasum * 1.0 / mrnacount
		sys.stderr.write("Total transcript length {}, average {:.2f} , longest is {}\n".format(mrnasum, avgmrna, longest_mrna) )

	# by default, should be the sum of the integers calculated for
	# either normal genes or when -g is used
	genesum = sum(genelengths)
	if genesum:
		sys.stderr.write("Total gene/transcript length (including overlaps) {}, longest is {}\n".format(genesum, max(genelengths) ) )
		nscaf = report_scaffold_stats(genesbyscaffold)
		geneOccupiedLength = get_occupied_length(genesbyscaffold)
		sys.stderr.write("Total non-redundant genic length  {}  {}\n".format(geneOccupiedLength, time.asctime() ) )

		if args.scaffolds: # scaffolds are given, meaning get gap information
			nscaf, scafn50, sclen, igaps, ggaps = find_gaps(args.scaffolds, exclusiondict, args.repeat, genesbyscaffold, args.use_original_scaffolds, args.verbose, args.stranded)

		if args.w: # special output for table
			# divide by one million for megabases
			if scafn50: # meaning if reassigned by find_gaps() args.scaffolds
				sys.stderr.write("{}\t{:.3f}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}\n".format(nscaf, scafn50/1000000.0, sclen/100000/10.0, exonOccupiedLength/100000/10.0, geneOccupiedLength/100000/10.0, igaps/100000/10.0, ggaps/100000/10.0) )
			else:
				sys.stderr.write("{}\t0\t0\t{:.1f}\t{:.1f}\t0\t0\n".format(nscaf, exonOccupiedLength/100000/10.0, geneOccupiedLength/100000/10.0 ) )
			if exonsum:
				sys.stderr.write("{}\t{}\t{}\t{}\n".format(nrexoncount, int(avgexon), introncount, int(avgintron) ) )
			else:
				sys.stderr.write("# WARNING: no exons found, cannot calculate exon and intron lengths\n")

		if args.print_intergenic: # treat intergenic regions as introns for most calculations
			intergenicsum, avgintergenic, intergeniccount = get_intron_stats( nametoscaffold, genesbyscaffold, args.print_intergenic, args.nested_genes, "intergenic intervals" )
		if args.repeat_masker:
			find_repeats( args.repeat_masker, genesbyscaffold, exclusiondict )

	# only for debugging
	if args.debug:
		# print some of the dictionaries
		for k,v in nametoscaffold.items():
			print("{}\t{}".format(k,v), file=sys.stdout)
		for k,v in exonboundaries.items():
			print("{}\t{}".format(k,v), file=sys.stdout)
		for k,v in genesbyscaffold.items():
			print("{}\t{}".format(k,v), file=sys.stdout)

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
