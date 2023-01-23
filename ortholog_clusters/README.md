# steps to make ortholog clusters #

### collection of datasets ###
A total of 35 candidate species were selected. In the end, 326130 total proteins were used for 15 species. Not all species were used, based on whether the proteomes were complete enough as measured by [BUSCO](https://busco.ezlab.org/) using the Metazoa set (954 genes, in table below), as well as the Eukaryotic set.

| species | complete | single | duplicated | fragmented | missing |
| :--- | --- | --- | --- | --- | --- |
| 01_COWC | 74.90% | 68.70% | 6.20% | 3.60% | 21.50% |
| 02_MBRE | 68.70% | 68.00% | 0.70% | 4.80% | 26.50% |
| 03_SROS | 66.10% | 65.50% | 0.60% | 6.90% | 27.00% |
| 04_AVAS | 80.30% | 78.60% | 1.70% | 2.30% | 17.40% |
| 05_AQUE | 89.90% | 84.40% | 5.50% | 4.20% | 5.90% |
| 06_EMUE | 70.50% | 58.10% | 12.40% | 6.90% | 22.60% |
| 07_SCIL | 72.50% | 70.40% | 2.10% | 10.80% | 16.70% |
| 08_HHON | 91.90% | 89.20% | 2.70% | 1.00% | 7.10% |
| 09_TADH | 91.70% | 90.70% | 1.00% | 2.70% | 5.60% |
| 10_BFLO | 88.70% | 88.10% | 0.60% | 5.10% | 6.20% |
| 11_HSAP | 100.00% | 98.80% | 1.20% | 0.00% | 0.00% |
| 12_CELE | 78.50% | 51.70% | 26.80% | 2.10% | 19.40% |
| DEM01_AASP | 9.10% | 7.20% | 1.90% | 10.30% | 80.60% |
| DEM02_FOSP | 61.80% | 52.90% | 8.90% | 21.00% | 17.20% |
| DEM03_LAIN | 65.80% | 46.50% | 19.30% | 12.40% | 21.80% |
| DEM04_LASP | 91.90% | 22.40% | 69.50% | 3.40% | 4.70% |
| DEM05_LISP | 87.30% | 32.40% | 54.90% | 6.90% | 5.80% |
| DEM06_MYSP | 48.70% | 44.20% | 4.50% | 25.80% | 25.50% |
| DEM07_STSP | 43.10% | 39.30% | 3.80% | 29.40% | 27.50% |
| DEM08_TESP | 41.40% | 36.40% | 5.00% | 23.90% | 34.70% |
| DEM09_TPAP | 46.00% | 41.00% | 5.00% | 26.50% | 27.50% |
| HEX01_AABY | 35.30% | 31.40% | 3.90% | 22.10% | 42.60% |
| HEX02_ASER | 76.20% | 48.50% | 27.70% | 7.20% | 16.60% |
| HEX03_BCYA | 70.20% | 57.20% | 13.00% | 9.90% | 19.90% |
| HEX04_CDIS | 73.50% | 52.10% | 21.40% | 7.40% | 19.10% |
| HEX05_CHSP | 71.20% | 43.20% | 28.00% | 10.20% | 18.60% |
| HEX06_EUSP | 79.20% | 67.80% | 11.40% | 3.70% | 17.10% |
| HEX07_FOCC | 79.00% | 69.30% | 9.70% | 4.10% | 16.90% |
| HEX08_FSIM | 17.90% | 15.20% | 2.70% | 22.00% | 60.10% |
| HEX09_LASP | 70.00% | 54.00% | 16.00% | 9.90% | 20.10% |
| HEX10_ROKI | 81.60% | 61.80% | 19.80% | 3.70% | 14.70% |
| HEX11_STET | 37.50% | 34.30% | 3.20% | 22.00% | 40.50% |
| HEX12_TRSP | 74.80% | 65.20% | 9.60% | 6.50% | 18.70% |
| HEX13_COSP | 25.60% | 22.70% | 2.90% | 21.90% | 52.50% |
| HEX14_WLEU | 65.40% | 36.40% | 29.00% | 14.00% | 20.60% |

### all-v-all BLAST and clustering ###
Using [diamond](https://github.com/bbuchfink/diamond) by [Buchfink et al 2017](https://www.nature.com/articles/nmeth.3176) instead of NCBI blastp, and using the script `makehomologs.py`, this requires `mcl` and some additional python libraries: [BioPython](https://biopython.org/) and [networkx](https://networkx.org/). Options refer to: 

* `-s 1` cluster must have at least 1 taxon
* `-z 2` cluster must have 2 or more
* `-L 40000` max protein length
* `-M 3000` cluster can have up to 3000 seqs (e.g. transposon clusters)
* `-H 1000` max seqs for a single taxon

```
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
```

This generates a log file `holozoa_allprots_v10.2022-03-21-165059.mh.log`, used to make a 3-page PDF output of several graphs `holozoa_allprots_v10.2022-03-21-165059.mh.pdf`, as well as other output files.

### alignments, tree inference, and domains ###
Align and make trees, using [MAFFT](https://mafft.cbrc.jp/alignment/server/) and [FastTree](www.microbesonline.org/fasttree/)

```
#!/bin/bash

NUMTH=4

# create MAFFT-einsi alignments
for file in *fasta ; do mafft --maxiterate 1000 --genafpair --thread $NUMTH $file >$file.aln ; done

# create FASTREE gene trees:
for file in *aln ; do fasttree -wag -gamma $file > $file.tre ; done
```

Get PFAM domains for all clusters using the [PFAM-A database](https://www.ebi.ac.uk/interpro/download/Pfam/) and [hmmscan](https://github.com/EddyRivasLab/hmmer):

```
#!/bin/bash

NUMTH=4
PFAM=/home/ubuntu/avas/ortholog_assignments/makehomologs/pfam_hmm/Pfam-A.hmm

# search PfamA with HMMSEARCH:
for file in *fasta ; do hmmscan --cpu $NUMTH --domtblout ${file%.fasta}.pfam.txt $PFAM $file ; done

# reformat output:
for file in *txt ; do tabularpfam.py $file > ${file%.txt}.tab ; done

mkdir original_output
mv *txt original_output
mkdir formatted_output
mv *tab formatted_output
```

### interactive viewing of clusters ###
Using [Rshiny](https://shiny.rstudio.com/), which can be run in Rstudio. This can, in theory, be used by anyone, or modified for your own clusters or genomes. Please contact [user WRF](https://github.com/wrf) with questions about implementation.

![ortholog_cluster_browser_app_screenshot_01.png](https://github.com/PalMuc/Aphrocallistes_vastus_genome/blob/main/ortholog_clusters/ortholog_cluster_browser_app_screenshot_01.png)

