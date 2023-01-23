# phylogeny of hexactinellids #
Most of the data is from the 2017 RV SONNE cruise as part of Project PoribacNewZ. The 454 transcriptome of [Euplectella aspergillum](https://marinegenomics.oist.jp/kairou/viewer/info?project_id=62) is from the [OIST Marine Genomics Unit](https://marinegenomics.oist.jp/gallery/gallery/index), and the three other transcriptomes come from [Whelan et al 2015](https://figshare.com/articles/dataset/Error_signal_and_the_placement_of_Ctenophora_sister_to_all_other_animals/1334306).

Most of the code can be found in [WRF's supermatrix repo](https://github.com/wrf/supermatrix).

Make overall supermatrix:

```
./collate_busco_results.py -s species_to_dir_list.tab -o combined_buscos/

for FILE in combined_buscos/*.faa ; do mafft $FILE > $FILE.aln ; done

~/git/supermatrix/join_alignments.py -a combined_buscos/*.aln -d "|" -u hexa_sponge_combined_buscos.aln
# Finished parsing 926 alignments for 22 total taxa
```

Sort and make plot for all taxa for all genes:

```
~/git/supermatrix/reorder_matrix_by_cov.py -a hexa_sponge_combined_buscos.aln -p hexa_sponge_combined_buscos.aln.partition.txt -o hexa_sponge_combined_buscos.sorted.aln
# alignment has 926 partitions with 526348 sites

~/git/supermatrix/check_supermatrix_alignments.py -a hexa_sponge_combined_buscos.sorted.aln -p hexa_sponge_combined_buscos.sorted.aln.partition.txt -m hexa_sponge_combined_buscos.occ_matrix.tab --percent
# Matrix has 101 (0.50%) complete and 9928 (48.73%) empty out of 20372 total genes
Rscript ~/git/supermatrix/draw_matrix_occupancy.R hexa_sponge_combined_buscos.occ_matrix.tab
```

Filter for 50% coverage and require certain low coverage taxa:

```
~/git/supermatrix/filter_supermatrix.py -a hexa_sponge_combined_buscos.aln -p hexa_sponge_combined_buscos.aln.partition.txt -r Hyalonema_populiferum Walteria_leuckarti Corbitella_sp Farrea_similaris -o hexa_sponge_combined_buscos.filt.aln
~/git/supermatrix/reorder_matrix_by_cov.py -a hexa_sponge_combined_buscos.filt.aln -p hexa_sponge_combined_buscos.filt.aln.partition.txt -o hexa_sponge_combined_buscos.filt.sorted.aln
~/git/supermatrix/check_supermatrix_alignments.py -a hexa_sponge_combined_buscos.filt.sorted.aln -p hexa_sponge_combined_buscos.filt.sorted.aln.partition.txt -m hexa_sponge_combined_buscos.filt.occ_matrix.tab --percent
Rscript ~/git/supermatrix/draw_matrix_occupancy.R hexa_sponge_combined_buscos.filt.occ_matrix.tab
```

Make draft tree:

`FastTree hexa_sponge_combined_buscos.filt.sorted.aln > hexa_sponge_combined_buscos.filt.sorted.tree`

