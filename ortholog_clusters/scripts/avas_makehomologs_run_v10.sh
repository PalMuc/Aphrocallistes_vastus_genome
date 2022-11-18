#!/bin/bash


NUMTH=40


SETE=holozoa_allprots_v10.fasta


# protein SETE:
diamond makedb --in $SETE --db $SETE
diamond blastp -q $SETE -d $SETE -p $NUMTH -e 0.01 --outfmt 6 > ${SETE%.fasta}.blastp.tab
makehomologs_g1.2.py -i ${SETE%.fasta}.blastp.tab -f $SETE -d _ -s 1 -z 2 -L 40000 -p 234 -M 3000 -H 1000 -T $NUMTH -o ${SETE%.fasta} -c 


tar czvf clusters_${SETE%.fasta}.tar.gz clusters_${SETE%.fasta}

mkdir holozoa_allprots_v10
mv *holozoa_allprots_v10* holozoa_allprots_v10/
