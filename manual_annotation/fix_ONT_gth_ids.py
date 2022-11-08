#!/usr/bin/env python

# fix_ONT_gth_ids.py

'''
fix_ONT_gth_ids.py BRAKER2_ONT-RNA_BRAKER2_ONT-RNA_protein_alignment_gth.gff3
'''

import sys

if len(sys.argv) < 2:
	sys.exit(__doc__)
else:
	for line in open(sys.argv[1]):
		line = line.strip()
		if line and line[0] != "#":
			lsplits = line.split('\t')
			scaffold = lsplits[0]
			scaf_id = scaffold.rsplit("_",1)[1]
			attributes = lsplits[8]
			new_attrs = attributes.replace("ID=","ID=s{}.".format(scaf_id)).replace("Parent=","Parent=s{}.".format(scaf_id))
			lsplits[8] = new_attrs
			newline = "{}\n".format( "\t".join(lsplits) )
			sys.stdout.write(newline)
