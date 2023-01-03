#!/usr/bin/env python

'''
fix_avas_rrna_gff.py Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_rRNA_revised.gff
'''

import sys
from collections import defaultdict

if len(sys.argv) < 2:
	sys.exit(__doc__)
else:
	rrna_gff = sys.argv[1]

	sys.stderr.write("# processing rRNA from {}\n".format(rrna_gff))
	linecounter = 0
	annot_counts = defaultdict(int)
	for line in open(sys.argv[1],'r'):
		if line.strip() and line[0]!="#":
			linecounter += 1
			lsplits = line.split("\t")
			attributes = lsplits[8].strip()
			feat_name = attributes.replace("Name=","").replace(" ",".")
			new_attr = "ID=Avas_rRNA_{}_{};{}".format(linecounter, feat_name, attributes)
			lsplits[2] = "rRNA"
			lsplits[8] = new_attr
			parent_line = "\t".join(lsplits)
			print(parent_line, file=sys.stdout)
			new_attr = "Parent=Avas_rRNA_{}_{}".format(linecounter, feat_name)
			lsplits[2] = "exon"
			lsplits[8] = new_attr
			exon_line = "\t".join(lsplits)
			print(exon_line, file=sys.stdout)

