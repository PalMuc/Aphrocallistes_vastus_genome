#!/usr/bin/env python
#
# v1.1 makehomologs.py
#
# code was forked from the transcriptome pipeline "Agalma" by the Dunn Lab
# https://bitbucket.org/caseywdunn/agalma/src/master/
# incorporates steps from the homologize and multalign pipeline stages
# with some modifications
#
# download mcl by:
# sudo apt install mcl
#
# and finally download diamond here (faster substitute for NCBI BLAST suite):
# https://github.com/bbuchfink/diamond
#
# for updates to networkx
# https://networkx.org/documentation/stable/reference/index.html

"""makehomologs.py v1.2  last modified 2022-03-17

  All files are written into either pwd or -o output_dir

  Starting from a fasta file of all sequences:
makehomologs.py -i allprots.fasta -o test -p 1234

  proteins are expected to be '|' delimited, in the format of
  Species_name|protein_ID

  otherwise delimiter can be changed with -d

Programs are:
1 - BLAST - diamond makedb + blastp;   expects input fasta
2 - filter - filter blast hits by networkx;  expects diamond/blast _hits.tab
3 - MCL - cluster genes by mcl;    expects allvall_edges_.abc
4 - refine - sort clusters;      expects mcl_gene_clusters_.abc
                                 returns 10 column file fasta_clusters._.tab
    filename  nseq  ntaxa  Nmin  Nmed  Nmax  Lmin  Lmean  Lmax  taxa-list

    or using option -c
    display counts of each taxon as a column (probably not useful beyond a dozen or so taxa)
    the tabular file will instead be named  fasta_clusters.H._.tab

  NOTE: diamond and mcl must be in $PATH

makehomologs.py -i lactobacillus_all_proteins.all_v_all.tab -f lactobacillus_all_proteins.fasta -p 234 -M 500 -T 4 -o orthologs_v1

  diamond step could be run as:

diamond-v2.0.13 blastp -q lactobacillus_all_proteins.fasta -d lactobacillus_all_proteins.fasta -o lactobacillus_all_proteins.all_v_all.tab

"""

import sys
import os
import traceback
import argparse
import subprocess
import time
import gzip
import multiprocessing
from collections import namedtuple, defaultdict
import networkx
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def safe_mkdir(path):
	if os.path.isfile(path):
		raise ValueError("{0} is a regular file, a directory is required".format(path))
	elif os.path.isdir(path):
		print("# The directory {0} already exists, and will be used".format(path), file=sys.stderr)
	else:
		print("# Making the directory {0}".format(path), file=sys.stderr)
		os.mkdir(path)


def logreport(loginput, wayout):
	# if input is string, then output like this:
	# This is a sample string for stderr Wed Oct 23 2013 16:31
	# or:
	# Wed Oct 23 2013 16:31
	# This is the sample string written to log
	if type(loginput) is str:
		print("{}\n{}".format(loginput, time.asctime()),  file=sys.stderr )
		print("#TIME-{}\n{}".format(time.asctime(), loginput), file=wayout )
	# otherwise if input is a list, then assume it is a list of arguments and print as such
	elif type(loginput) is list:
		print("Making system call:\n{}".format(" ".join(loginput)),  file=sys.stderr )
		print("#TIME-{}\n#CMD  {}".format(time.asctime(), " ".join(loginput)), file=wayout )


def prepare_blast(output_dir, fasta, wayout):
	"""make diamond database for blastp, and run all against all, then return the output filename"""
	print("Prepare all-by-all DIAMOND database", file=sys.stderr )

	BlastDBArgs = ["diamond makedb", "--in", fasta, "--db", fasta]
	logreport(BlastDBArgs,wayout)
	subprocess.call(BlastDBArgs)

	diamondresult = "{}.diamond_all_v_all.tab".format(output_dir)
	diamondBlastpArgs = ["diamond blastp", "-q", fasta, "-d", fasta, "-o", diamondresult]
	logreport(diamondBlastpArgs,wayout)
	subprocess.call(diamondBlastpArgs)

	return diamondresult


def get_seq_lengths(sequences):
	'''from a fasta file, return a dictionary where protein ID is the key and length is the value'''
	seqlendict = {}
	print( "# Parsing sequences from {}  {}".format(sequences, time.asctime() ), file=sys.stderr)
	for seqrec in SeqIO.parse(sequences,'fasta'):
		seqlendict[seqrec.id] = len(seqrec.seq)
	return seqlendict


def parse_edges(output_dir, blast_hits, min_overlap, min_bitscore, min_cluster_size, species_delimiter, verbose, wayout, seqfile=None):
	"""Parse BLAST hits into edges weighted by bitscore"""

	## step 2 ##
	############ ############
	# expects input as blast_hits.tab outfmt 6 from blast or diamond output

	# should deal with:
	# full length homologs - length and bitscore required
	# partial true homologs - must both hit non-overlapping parts of the same protein
	# splice variants
	# multiple paralogs - should be equally distant from homolog, thus take both

	nGr = networkx.Graph()
	edge_file = "allvall_edges_{}.abc".format(output_dir)

	seqlendict = get_seq_lengths(seqfile)
	protein_counts_by_sp = defaultdict(int) # key is species code, value is count
	for seqid in seqlendict.keys():
		taxon_identifier = seqid.split(species_delimiter)[0]
		protein_counts_by_sp[taxon_identifier] += 1

	logreport( "# Counted {} total proteins for {} species  {}".format( len(seqlendict), len(protein_counts_by_sp), time.asctime() ), wayout )
	logreport( "\n".join( ["T_genes\t{}\t{}".format(k, v) for k,v in sorted(protein_counts_by_sp.items(),key=lambda x:x[0]) ] ) , wayout )

	#max overlap is inverse of min_overlap
	max_overlap = float(1/min_overlap)

	seqs_w_blast_hit = {} # key is seq ID, value is True
	edge_type_counter = {'all': 0, 'non-self': 0, 'passed-min-overlap': 0, 'passed-max-overlap':0, 'passed-bitscore': 0, 'retained-edges': 0}


	if blast_hits.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		_opentype = gzip.open
		print("# Starting parsing on {} as gzipped".format( blast_hits ), file=sys.stderr)
	else: # otherwise assume normal open for fasta format
		_opentype = open
		print("# Starting parsing on {}".format( blast_hits ), file=sys.stderr)
	with _opentype(blast_hits,'rt') as bh:
		for line in bh:
			# fields:
			blastsplits = line.strip().split('\t')
			id_from = blastsplits[0]
			id_to = blastsplits[1]
			bitscore = blastsplits[11]
			qlen = seqlendict[id_from]
			slen = seqlendict[id_to]


			bitscore = float(bitscore)
			length = float(slen) / float(qlen)
			# Correct for nucleotide vs. amino acid length
			# only matters if program was tblastx
			#if seq_type == 'nucleotide':
			#	slen *= 3.0
			# Filter out self hits, low scoring hits, and short hits
			edge_type_counter['all'] += 1
			# discard self hits, meaning seq hits itself, probably with maximum bitscore
			if id_from == id_to:
				selfbitscore = bitscore
			else:
				edge_type_counter['non-self'] += 1
				# discard hits shorter than overlap cutoff
				if length >= min_overlap:
					edge_type_counter['passed-min-overlap'] += 1
					if length <= max_overlap:
						edge_type_counter['passed-max-overlap'] += 1
						# discard hits where bitscore / query length is too low
						# bits/subject length tends to find too many short sequences, which cannot be assigned efficiently
						# so instead query length is used, as this can only tolerate long and good hits
						bpl = bitscore/float(qlen)
						if bpl >= min_bitscore:
							edge_type_counter['passed-bitscore'] += 1
							# If an edge already exists between the nodes, update its score to be the max
							# this will occur for hits back, use the highest score for edge weight
							# this may cause problems for large genes where multiple hsps are found 2014-08-16
						#	if nGr.has_edge(id_from, id_to):
						#		e = nGr.edge[id_from][id_to]
						#		e['score'] = e.get('score')+bitscore
							# otherwise create new node
						#	else:
								#dg.add_node(id_from)
								#dg.add_node(id_to)
							nGr.add_edge(id_from, id_to, score=bitscore)
							seqs_w_blast_hit[id_from] = True
							seqs_w_blast_hit[id_to] = True
						elif verbose:
							print( "BITS %s %s %.2f" % (id_from, id_to, bpl) , file=sys.stderr)
					elif verbose:
						print( "LONG %s %s %.2f" % (id_from, id_to, length) , file=sys.stderr)
				elif verbose:
					print( "SHORT %s %s %.2f" % (id_from, id_to, length) , file=sys.stderr)

	blast_hits_by_sp = defaultdict(int)
	for seqid in seqs_w_blast_hit.keys():
		taxon_identifier = seqid.split(species_delimiter)[0]
		blast_hits_by_sp[taxon_identifier] += 1
	logreport( "\n".join( ["N_bhits\t{}\t{}".format(k, v) for k,v in sorted(blast_hits_by_sp.items(),key=lambda x:x[0]) ] ) , wayout )

	# the whole purpose of the subgraph was to allow for counting of nodes
	# however the connected_components function gives a list of lists
	# which can easily be added to a list and removed after iteration
	nnodes = {'nodes': 0, 'filtered-nodes': 0}
	logreport( "# Counting nodes in subgraphs" , wayout )
	nodestoprune = []
	for cc in networkx.connected_components(nGr):
		nodecount = len(cc)
		nnodes['nodes'] += nodecount
		# remove nodes in connected components if the number of nodes is fewer than specified value
		# this has no effect on very large clusters
		if nodecount < min_cluster_size: # default is < 4
			nodestoprune.extend(cc)
			# this will count the contribution of each size below the min
			nnodes[str(nodecount)] = nnodes.get(str(nodecount),0)+1
	print( "# Writing final graph file {}".format( edge_file ) , file=sys.stderr)
	logreport( "# Found {} subgraphs, pruning small subgraphs ({} nodes)".format( nnodes.get("nodes"), len(nodestoprune) ) , wayout )
	nGr.remove_nodes_from(nodestoprune)

	total_node_count = 0
	# after nodes are removed, write remaining edges to file
	print( "# Writing final graph file {}".format( edge_file ) , file=sys.stderr)
	with open(edge_file, 'w') as f:
		nnodes['filtered-nodes'] = nGr.number_of_nodes()
		nodecount = nGr.number_of_nodes()
		total_node_count += nodecount
		print( "# Writing edges to {}".format( edge_file ) , file=sys.stderr)
		for id_from, id_to in networkx.edges(nGr, nbunch=None):
			edge_type_counter['retained-edges'] += 1
			# syntax changed in networkx
			#print( "{}\t{}\t{}".format( id_from, id_to, nGr.edge[id_from][id_to]['score']), file=f )
			print( "{}\t{}\t{}".format( id_from, id_to, nGr[id_from][id_to]['score']), file=f )

	# report final counts to log
	logreport( "\n".join( [ "F_edge\t{}\t{}".format(k,v) for k,v in edge_type_counter.items() ] ) , wayout )

	# report final node counts to log, including number of filtered nodes by size
	logreport( "\n".join( ["N_nodes\t{}\t{}".format(k,v) for k,v in nnodes.items() ] ) , wayout )
	if not total_node_count:
		print("# ERROR: no sequences were written to the FASTA file", file=sys.stderr)
	return edge_file



def mcl_cluster(edge_file, cluster_inflation, outdir, threadcount, wayout):
	"""run mcl and return file name of the output"""

	### step 3 ###
	# assumes input is allvall_edges_.abc

	print( "Run mcl on all-by-all graph to form gene clusters", file=sys.stderr)
	cluster_file = "mcl_gene_clusters_{}.abc".format(outdir)
	mlcArgs = ['mcl', edge_file, "--abc", "-I", cluster_inflation, '-o', cluster_file, "-te", threadcount]
	logreport(mlcArgs,wayout)
	subprocess.call(mlcArgs)
	return os.path.abspath(cluster_file)


def refine_clusters(outdir, clusterfile, min_clust_size, min_taxa, max_sequences, max_length, fastadb, max_homologs, species_delimiter, get_sp_counts, verbose, wayout):
	print( "Select a cluster for each homologize component that meets size, sequence length, and composition requirements", file=sys.stderr)

	##########
	# step 4 #
	##########
	# assumes input is mcl_gene_clusters_.abc

	cluster_dir = "clusters_{}".format(outdir)
	safe_mkdir(cluster_dir)

	clustsizehist = defaultdict(int) # key is number of sequences in cluster, value is count
	taxasizehist = defaultdict(int) # key is number of unique taxa per cluster, value is count
	seqs_type_counter = {'total-clusters': 0, 'seq_under_max_len': 0, 'max_sequences': 0, 'min_taxa': 0, 'taxa_mean': 0, 'single_copy_gene':0 , "paralogous_gene":0, "single_species":0 }
	seqs_type_descs = {'total-clusters': "total number of clusters" , 'seq_under_max_len': "total sequences shorter than length cutoff {}".format(max_length) , 'max_sequences': "clusters with fewer than {} sequences".format(max_sequences) , 'min_taxa': "clusters with more than {} taxa".format(min_taxa) , 'taxa_mean': "clusters with fewer than {} seqs per taxon".format(max_homologs), 'single_copy_gene':"clusters of a single copy ortholog" , "paralogous_gene":"clusters with more than one copy in any taxon" , "single_species":"cluster has only one taxon" }

	seqcounts = 0
	components = defaultdict(list) # key is cluster ID, value is list of seq IDs for that cluster

	# read fasta proteins, report total count
	dbdictionary = SeqIO.to_dict(SeqIO.parse(fastadb,"fasta"))
	logreport( "# fasta file contains {} records".format( len(dbdictionary) ) , wayout )
	# get counts of each species, mostly for keeping track of all species in the working dataset
	species_name_dict = defaultdict(int) # key is species, value is count
	for seqrec in dbdictionary:
		taxon_identifier = dbdictionary[seqrec].id.split(species_delimiter)[0]
		species_name_dict[taxon_identifier] += 1
	species_list = sorted(species_name_dict.keys()) # list of all species by their name/codes

	# read in cluster file and generate components contain each seq in the cluster
	with open(clusterfile, 'rt') as cl:
		# cluster counts begin at 1, not 0, for convenience of sorting later
		cluster_id = 1
		for line in cl:
			seqs_type_counter['total-clusters'] += 1
			cluster = line.rstrip().split()
			csize = len(cluster)
			clustsizehist[csize] += 1
			if csize >= min_clust_size:
				seqcounts += csize
				for seq_id in cluster:
					if len(dbdictionary[seq_id].seq) <= max_length:
						seqs_type_counter['seq_under_max_len'] += 1
						components[cluster_id].append(SeqRecord(dbdictionary[seq_id],id=seq_id))
				cluster_id += 1
	print( "# cluster size stats (including those removed by -z), as N sequences:number of clusters", file=sys.stderr)
	logreport( "\n".join( ["N_seqs\t{}\t{}".format(k,v) for k,v in sorted(clustsizehist.items(), key=lambda x: x[0], reverse=True) ] ) , wayout )

	occurrence_of_taxa = defaultdict(int) # key is taxon, value is count
	ortho_taxa = defaultdict(int) # key is taxon, value is count but only for 1-1 orthologs

	# count the number of components to get the number of zeroes
	xzeroes = str(len(str(len(components))))

	# write all cluster filenames to a log, for the next step of alignment
	if get_sp_counts:
		fastaclusters_list = "fasta_clusters.H.{}.tab".format(outdir)
		fastaclusters_headings = ["filename", "num_sequences", "num_taxa", "min_per_taxon", "med_per_taxon", "max_per_taxon", "min_length", "mean_length", "max_length"]+ species_list
	else:
		fastaclusters_list = "fasta_clusters.{}.tab".format(outdir)
		fastaclusters_headings = ["filename", "num_sequences", "num_taxa", "min_per_taxon", "med_per_taxon", "max_per_taxon", "min_length", "mean_length", "max_length", "taxaList" ]
	print( "# writing cluster information to {} \n# headings are:\n{}".format( fastaclusters_list, "\t".join(fastaclusters_headings) ), file=sys.stderr)
	# open tabular output of cluster names, and stats
	with open(fastaclusters_list, 'w') as fc:
		# if listing species with counts, add header line, which may disrupt some downstream analysis that does not expect headers
		if get_sp_counts:
			print("\t".join(fastaclusters_headings) , file=fc)
		# Loop over the components and write the fasta files
		for comp_id, records in components.items():
			# must reset the taxa_count dictionary for each component
			taxa_count = defaultdict(int)
			# Apply filter on cluster size, which is not done before as long seqs can be filtered
			csize = len(records) # number of sequences in cluster
			seq_length_list = [] # length in AAs of each seq in cluster

			if csize <= max_sequences:
				seqs_type_counter['max_sequences'] += 1

				# Count the number of times a taxon is repeated in the cluster
				for record in records:
					seq_length_list.append( len(record.seq) )
					taxon_identifier = record.id.split(species_delimiter)[0] # default split is at first |
					occurrence_of_taxa[taxon_identifier] += 1
					taxa_count[taxon_identifier] += 1
					#taxa_count[record.id[0:2]] = taxa_count.get(record.id[0:2], 0) + 1
				# Compute the mean number of repetitions
				count_by_taxa = taxa_count.values()
				ntaxa = len(count_by_taxa)
				taxa_mean = sum(count_by_taxa) / float(ntaxa)
				taxasizehist[ntaxa] = taxasizehist.get(ntaxa,0) + 1
				if ntaxa == 1: # cluster has only single species
					seqs_type_counter["single_species"] += 1
				# Require at least min_taxa, and that the mean reptitions is less than the value of max_homologs.
				if ntaxa >= min_taxa:
					seqs_type_counter['min_taxa'] += 1
					if taxa_mean < max_homologs:
						seqs_type_counter['taxa_mean'] += 1

						# define each cluster as orthologs of 1-1 or paralogs of any taxa
						if csize == ntaxa:
							seqs_type_counter['single_copy_gene'] += 1
							clustoutname = 'homologs_{1}_{2:0{0}}.fasta'.format(xzeroes, outdir, comp_id)
							for taxon in taxa_count.keys():
								ortho_taxa[taxon] += 1
						else:
							seqs_type_counter['paralogous_gene'] += 1
							clustoutname = 'paralogs_{1}_{2:0{0}}.fasta'.format(xzeroes, outdir, comp_id)

						# write fasta format of each cluster
						clust_fasta_file = os.path.join(cluster_dir, clustoutname)
						with open(clust_fasta_file, 'w') as fo:
							for record in records:
								fo.write("%s" % dbdictionary[record.id].format("fasta"))

						# write cluster stats to the big table
						if get_sp_counts:
							cl = [clust_fasta_file, csize, ntaxa, min(count_by_taxa), sorted(count_by_taxa)[(ntaxa)//2], max(count_by_taxa), min(seq_length_list), sum(seq_length_list)/len(seq_length_list), max(seq_length_list), "\t".join( [ str(taxa_count.get(sp,0)) for sp in species_list] ) ]
						else:
							cl = [clust_fasta_file, csize, ntaxa, min(count_by_taxa), sorted(count_by_taxa)[(ntaxa)//2], max(count_by_taxa), min(seq_length_list), sum(seq_length_list)/len(seq_length_list), max(seq_length_list), ",".join(list(taxa_count.keys())) ]
						# should be 10+ columns of:
						# file  nseq ntaxa min med max min mean max  taxa-list
						print( "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1f}\t{}\t{}".format( *cl ), file=fc )

					else: # too many paralogs in some/many species, meaning taxa_mean < max_homologs
						logreport( "PARALOGS-{}\t{:.2f}\tremoved".format(comp_id, taxa_mean) , wayout )
						#print( "PARALOGS-{}\t{:.2f}\tremoved".format(comp_id, taxa_mean), file=sys.stderr)
				else: # not enough taxa
					pass
					#logreport( "SMALL-{}\t{}\tremoved".format(comp_id, ntaxa) , wayout ) # prints thousands of lines
			else: # too many sequences, meaning csize > max_sequences
				logreport( "LARGE-{}\t{}\tremoved".format(comp_id, csize) , wayout )
				#print( "LARGE-{}\t{}\tremoved".format(comp_id, csize), file=sys.stderr)

	# report final counts to log
	print( "# cluster stats:", file=sys.stderr)
	logreport( "\n".join( ["F_clust\t{}: {}\t{}".format(k, seqs_type_descs.get(k), v) for k,v in seqs_type_counter.items() ] ) , wayout )

	logreport("# found {} total taxa".format(len(occurrence_of_taxa)) , wayout )
	logreport( "\n".join( ["N_genes\t{}\t{}".format(k, v) for k,v in sorted(occurrence_of_taxa.items(), key=lambda x:x[0]) ] ) , wayout )
	logreport( "\n".join( ["N_ortho\t{}\t{}".format(k, v) for k,v in sorted(ortho_taxa.items(), key=lambda x:x[0]) ] ) , wayout )

	# report dictionary histogram as string
	print( "# taxon count stats, as N taxa:number of clusters with N", file=sys.stderr)
	logreport( "\n".join( ["N_taxa\t{}\t{}".format(k,v) for k,v in sorted(taxasizehist.items(), key=lambda x: x[0], reverse=True) ] ) , wayout )

	return fastaclusters_list


def get_threads(threads):
	if threads:
		threadcount = threads
	else:
		print( "Detecting processors...", file=sys.stderr)
		threadcount = str(multiprocessing.cpu_count())
	print( "Using {} threads".format(threadcount), file=sys.stderr)
	return threadcount


def main(argv):
	if not len(argv):
		argv.append("-h")

	program_list = "1234"

	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-i','--input', help="set of fasta sequences, or blast hits, or edges", required=True)
	parser.add_argument('-f','--fasta', help="original fasta sequences, if starting after blast")
	parser.add_argument('-o','--output-dir', help="directory of output files", default='', required=True)
	parser.add_argument('-p','--program', help="string for programs used: '{}'".format(program_list), default="234")

	parser.add_argument('-z','--min_cluster_size', type=int, default=4, help="""
	Minimum number of sequences (nodes) to be retained as a homolgous gene cluster. [4]""")
	parser.add_argument('-s','--min_taxon_count', type=int, default=2, help="""
	Minimum number of different species to keep the cluster. [2]""")

	parser.add_argument('-b','--min_bitscore', type=float, default=0.7, help="""
	Remove BLAST hit edges with a bitscore/length ratio below min_bitscore. [0.7]""")
	parser.add_argument('-m','--min_overlap', type=float, default=0.5, help="""
	Filter out BLAST hit edges with a max HSP length less than this fraction
	of target length. [0.5]""")

	parser.add_argument('-I','--cluster-inflation', default="2.1", help="""
	MCL cluster inflation parameter -I, ranging from appx 1-5, 
	higher numbers make more splits/clusters. [2.1]""")

	parser.add_argument('-M','--max_sequences', type=int, default=300, help="""
	Maximum number of sequences to allow in a cluster. [300]""")
	parser.add_argument('-L','--max_length', type=int, default=15000, help="""
	Maximum length for sequences in a cluster [15000].""")
	parser.add_argument('-H','--max_homologs', type=int, default=6, help="""
	Maximum average number of sequences from one taxa. [6]""")
	parser.add_argument('-d','--species-delimiter', default='|', help="delimiter of sequence names in clustering step, of species|sequence_name. default [|]")

	parser.add_argument('-T','--threads', metavar='N', help="number of threads/processors, default: autodetect all")
	parser.add_argument('-c','--count-taxa', action="store_true", help="print counts of each taxon for all clusters in fasta_clusters.output file")
	parser.add_argument('-l','--log-dir', help="Optional directory of log files", default='./')
	parser.add_argument('-v','--verbose', action="store_true", help="display more output to stderr")
	args = parser.parse_args(argv)

	print("# Starting process:  {}".format( time.asctime() ), file=sys.stderr)
	startclock=time.asctime()
	starttime= time.time()

	## must change num_threads if running on a different system
	threadcount = get_threads(args.threads)
	
	# parse program string to get operations
	program_string = args.program
	get_programs = [bool(args.program.count(x)) for x in program_list]
	do_blast, do_parse, do_cluster, do_refine = get_programs

	# filename processing for output log
	log_file = "{0}.{1}.mh.log".format(args.output_dir, time.strftime("%Y-%m-%d-%H%M%S"))
	log_path = os.path.join(args.log_dir, log_file)
	print("# Log files being written to {}:  {}".format(log_path , time.asctime() ), file=sys.stderr)
	wayout = open(log_path, 'w')

	# write all the arguments from the command line to file
	logreport(" ".join(sys.argv), wayout )

	# get version info and write to the log
	minorversion = sys.version_info[1]
	logreport("# Using Python version: %s" % (".".join([str(i) for i in sys.version_info[0:3]])), wayout )

	# standard output of file counts
	logreport("# {} input file detected".format(args.input), wayout )
	if args.fasta:
		logreport("# {} used as fasta database".format(args.fasta), wayout )

	if args.min_taxon_count < 1:
		print("# WARNING: min taxon count -s {} cannot be less than 1, forcing to 1".format( args.min_taxon_count ), file=sys.stderr)
		args.min_taxon_count = 1
	if args.min_cluster_size < 1:
		print("# WARNING: min cluster size -z {} cannot be less than 1, forcing to 1".format( args.min_cluster_size ), file=sys.stderr)
		args.min_cluster_size = 1
	if args.min_taxon_count > args.min_cluster_size:
		print("# WARNING: min cluster size -z {1} cannot be less than min taxon count -s {0}, forcing cluster size to {0}".format( args.min_taxon_count, args.min_cluster_size ), file=sys.stderr)
		args.min_cluster_size = int(args.min_taxon_count)

	# beginning of main processes
	try:
		if do_blast:
			blast_hits = prepare_blast( args.output_dir, args.input, wayout )
			fasta_db = args.input
		else:
			blast_hits = args.input
			fasta_db = args.fasta
		if do_parse:
			edge_file = parse_edges(args.output_dir, blast_hits, args.min_overlap, args.min_bitscore, args.min_cluster_size, args.species_delimiter, args.verbose, wayout, args.fasta)
		else:
			edge_file = args.input
		if do_cluster:
			cluster_file = mcl_cluster(edge_file, args.cluster_inflation, args.output_dir, threadcount, wayout)
		else:
			cluster_file = args.input
		if do_refine:
			fastafiles = refine_clusters(args.output_dir, cluster_file, args.min_cluster_size, args.min_taxon_count, args.max_sequences, args.max_length, fasta_db, args.max_homologs, args.species_delimiter, args.count_taxa, args.verbose, wayout)
		#else:
		#	fastafiles = args.input

	except:
		#prints the error to the log file, as of 18-mar-2013
		traceback.print_exc(file=sys.stderr)
		traceback.print_exc(file=wayout)

	# final output
	logreport("#TIME-Process completed in %.1f minutes" % ((time.time()-starttime)/60), wayout )
	print( "Done:  {}".format(time.asctime()) , file=sys.stderr)
	print( "## Great success!!:  {}".format(time.asctime()) , file=wayout)
	wayout.close()
	return 1

if __name__ == "__main__":
	main(sys.argv[1:])
