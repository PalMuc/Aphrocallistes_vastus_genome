#!/usr/bin/env python
#
# get_read_skip_from_bam.py
# original version made for Schultz et al 2021 G3
# 
# current version available at
# https://bitbucket.org/wrf/sequences/src/master/

'''
get_read_skip_from_bam.py  last modified 2020-09-03
    from a SAM or BAM file
    extract length of skipped bases at beginning or end of a long read

~/samtools-1.9/samtools view UCSC_Hcal_v1_B1_LR.sorted.bam | get_read_skip_from_bam.py - > UCSC_Hcal_v1_B1_LR.sorted.read_skip.tab

    creates a 4-column table, of line number, starting -S, ending -S, and sequence of -S

cut -f 4 UCSC_Hcal_v1_B1_LR.sorted.read_skip.tab | sort | uniq -c | sort -nr

1845104 0
 354257 CCCTCAAAGTTTGAAAAGTTGTGATGAAATTTGTTTAATTAAAC
 321762 GGGAGTTTCAAACTTTTCAACACTACTTTAAACAAATTAATTTG

'''

import sys
import argparse
import re
import time

def cigar_to_list(cigarstring):
	'''from a cigar string, return a list of each number letter pair'''
	cigar_list = [] # build growing list with each RE
	for rematch in re.finditer("\d+",cigarstring):
		rematch_start = rematch.start()
		rematch_end = rematch.end() + 1
		number_letter_pair = cigarstring[rematch_start:rematch_end]
		cigar_list.append(number_letter_pair)
	return cigar_list

def use_skipped_sequence(motif):
	'''check if motif is polyA or polyT, and return boolean'''
	motif_len = len(motif)
	T_percent = motif.count("T") * 1.0 / motif_len
	A_percent = motif.count("A") * 1.0 / motif_len
	if T_percent >= 0.9 or A_percent >= 0.9:
		return False
	else:
		return True

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('input_file', type = argparse.FileType('rU'), default = '-', help="SAM file, - for stdin from BAM")
	parser.add_argument('-n', '--target-min', type=int, default=42, help="print out leader sequences with S value of at least N [42]", metavar="N")
	parser.add_argument('-N', '--target-max', type=int, default=45, help="print out leader sequences with S value of at most N [45]", metavar="N")

	args = parser.parse_args(argv)

	linecounter = 0
	cigar_counter = 0

	match_target_count = 0
	polyA_count = 0
	no_cigar_count = 0

	sys.stderr.write("# Reading {}, tracking S of {}-{}bp  {}\n".format(args.input_file.name, args.target_min, args.target_max, time.asctime() ) )
	for line in args.input_file:
		linecounter += 1
		lsplits = line.split("\t")
		if len(lsplits) < 12:
			continue
		cigar_string = lsplits[5]
		cigar_list = cigar_to_list(cigar_string)
		if len(cigar_list)==0:
			no_cigar_count += 1
			continue
		five_pr_skip = cigar_list[0][0:-1] if cigar_list[0][-1]=="S" else 0
		three_pr_skip = cigar_list[-1][0:-1] if cigar_list[-1][-1]=="S" else 0
		cigar_counter += 1

		leader_seq = "0"

		seq_string = lsplits[9]
		if args.target_min <= int(five_pr_skip) <= args.target_max:
			match_target_count += 1
			skipped_seq = seq_string[0:int(five_pr_skip)]
			if use_skipped_sequence(skipped_seq):
				leader_seq = skipped_seq
			else:
				polyA_count += 1
		if args.target_min <= int(three_pr_skip) <= args.target_max:
			match_target_count += 1
			seq_string_len = len(seq_string)
			skipped_seq = seq_string[ (seq_string_len - int(three_pr_skip)):]
			if use_skipped_sequence(skipped_seq):
				leader_seq = skipped_seq
			else:
				polyA_count += 1

		sys.stdout.write("{}\t{}\t{}\t{}\n".format( linecounter, five_pr_skip, three_pr_skip, leader_seq ) )

	sys.stderr.write("# Counted {} lines for {} reads  {}\n".format(linecounter, cigar_counter, time.asctime() ) )
	if no_cigar_count:
		sys.stderr.write("# {} had no CIGAR string\n".format( no_cigar_count ) )
	if match_target_count:
		sys.stderr.write("# {} reads had either end matched to target leader length\n".format(match_target_count) )
	if polyA_count:
		sys.stderr.write("# {} were likely polyA tails, and were ignored\n".format( polyA_count ) )

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
