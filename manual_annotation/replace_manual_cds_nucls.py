#!/usr/bin/env python

"""
./replace_manual_cds_nucls.py Avas.v1.28_annotations.gff manual_genes_to_map.fasta Avas.v1.28_annotations.nucl.fasta > Avas.v1.28_annotations.corr_nucl.fasta
"""

import sys
import re
from Bio import SeqIO

if len(sys.argv)<2:
	sys.exit(__doc__)
else:
	tx_to_source_id = {}
	rRNA_ids = {} # remove CDS of all of these
	annot_file = sys.argv[1]

	sys.stderr.write("# reading {} for seqs to remove\n".format(annot_file) )
	for line in open(annot_file,'r'):
		line = line.strip()
		if line and line[0]!='#':
			lsplits = line.split('\t')
			feature = lsplits[2]
			attributes = lsplits[8]
			try:
				gff_id = re.search("ID=([^;\s]+);", attributes).group(1)
			except AttributeError:
				continue
			if feature=="rRNA":
				rRNA_ids[gff_id] = True
			elif feature=="transcript":
				source_id = re.search("source_ID=([^;\s]+)", attributes).group(1)
				tx_to_source_id[gff_id] = source_id
	sys.stderr.write("# found {} rRNA IDs\n{}\n".format( len(rRNA_ids),rRNA_ids.keys() ) )
	sys.stderr.write("# found {} source IDs\n".format( len(tx_to_source_id) ) )

	manual_cds_file = sys.argv[2]
	sys.stderr.write("# reading {} for seqs to replace\n".format(manual_cds_file) )
	manual_gene_dict = SeqIO.to_dict(SeqIO.parse(manual_cds_file,"fasta"))
	sys.stderr.write("# found {} for seqs to replace\n".format(len(manual_gene_dict)) )

	nucl_file = sys.argv[3]
	for seqrec in SeqIO.parse(nucl_file,"fasta"):
		seqid = str(seqrec.id)
		if seqid in rRNA_ids:
			print( "# removing rRNA seq {}".format(seqid), file=sys.stderr)
			continue
		source_id = tx_to_source_id.get(seqid,None)
		if source_id in manual_gene_dict:
			print( "# replacing {} with {}".format(seqid, source_id), file=sys.stderr)
			seqrec.seq = manual_gene_dict.get(source_id).seq
		sys.stdout.write(seqrec.format("fasta"))
	sys.stderr.write("# done\n" )



