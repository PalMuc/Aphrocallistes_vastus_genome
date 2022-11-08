#!/usr/bin/env python

'''
./prepare_names_for_spreadsheet.py BRAKER2_ONT-RNA_minimap2_augustus.hints.gff3 > names_for_spreadsheet
'''


import sys
from collections import defaultdict

scaf_to_genes = defaultdict(list)
gc = 0

for line in open(sys.argv[1],'r'):
	lsplits = line.strip().split("\t")
	scaf = lsplits[0]
	attributes = lsplits[8]
	feature = lsplits[2]
	if feature=="gene":
		gc += 1
		short_scaf = scaf.replace("Aphrocalllistes_vastus_HiC-","")
		gene_id = attributes.split(";")[0].replace("ID=","")

		scaf_to_genes[short_scaf].append(gene_id)
print >> sys.stderr, "total is {}".format(gc)

csum = 0

for scaf in sorted(scaf_to_genes.keys()):
	csum += len(scaf_to_genes.get(scaf))
	print >> sys.stderr, "{} {} {}".format(scaf, len(scaf_to_genes.get(scaf)), csum)
#	for gene in sorted(scaf_to_genes.get(scaf), key=lambda x: int(x.split("gene")[1])):
	for gene in sorted(scaf_to_genes.get(scaf), key=lambda x: int(x.split("g")[1])):
		print >> sys.stdout, "{}\t{}".format(scaf, gene)
