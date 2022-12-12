# steps to count nested genes #
Nested genes refer to genes that are inside the intron of another gene, and could be on the same or opposite strand. These can be spliced or not, but are typically single exon.

```
 XXXXXXX>---------XXXXX>--------------XXXXXXX>  host gene, forward strand
           XXX>                                nested gene, same strand
                          <XXX---XXX           nested gene, opposite strand
```

### count nested with each genome ###
Nested gene count and base count were made using a strategy similar to that used by [Schultz 2021](https://doi.org/10.1093/g3journal/jkab302), with a few differences:

* [Schultz's code](https://github.com/conchoecia/chep/blob/master/scripts/gff_to_intron_bed.py) filtered to remove the largest introns, potentially due to mismapping)
* They also did the calcuations of transcript overlap by removing 5% of the length of both ends, again potentially reducing possibilities of original mismapping in the gene models, but still allows a small overlap
* Here, our script `gtfstats.py` required that a gene interval is completely within an intron interval from any strand or gene on that scaffold (i.e. that multiple isoforms may have similar introns, and gene could fit within any of them)

```
gtfstats.py -i v1.29/Avas.v1.29_annotations.fr.gff -s v1.29/Avas.v1.29_scaffolds.fasta -G -N -M Parent
gtfstats.py -N -i JAKMXF01.1.gbff.gff -s JAKMXF01.1.fsa_nt.fasta -w

gtfstats.py -i Emu_augustus_vs_nb_sysnames.gff --exon-counter -N
gtfstats.py -i Amphimedon_queenslandica.Aqu1.45.cleaned.gff -N -s Amphimedon_queenslandica.Aqu1.dna.toplevel.fa -w -G -M Parent
gtfstats.py -i sycon.cds.gff3 -N -G -M Parent

gtfstats.py -i ML2.2.gff3 -N -s MlScaffold09.nt.fa -G -M Parent
gtfstats.py -i Hcv1av93_release/Hcv1av93.gff -s UCSC_Hcal_v1.fa -w -N -G -M Parent

gtfstats.py -N -O -c -i Trichoplax_scaffolds_JGI_AUGUSTUS_no_comment.gff
gtfstats.py -N -i tracks/PlacoH13_BRAKER1_augustus_no_comment.gff -c
```

Most of the differences are due to annotation method. This is pronounced as [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus) contains the option:

```
--singlestrand=true
  predict genes independently on each strand, allow overlapping genes on opposite strands
  This option is turned off by default.
```

### Table 1: nested genes across non-bilaterian genomes ###
| species | N scaffolds | total Mbp | N genes | genic Mbp | exon Mbp | N nested | nested bp | nested exon | annot version | asm version | annot method |
| :--- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- |
| Aphrocallistes vastus | 186 | 80.320 | 20060 | 49.161 | 27.127 | 2432 | 3105825 | 2083470 | v1.29 | v1.29 | manual |
| Oopsacas minuta | 365 | 61.460 | 16340 | 29.089 | 20.650 | 83 | 60891 | 44676 | v1 | v1 | BRAKER |
| Ephydatia muelleri | 1444 | 322.620 | 39329 | 179.121 | 51.532 | 0 | 0 | 0 | v1 | v1 | AUGUSTUS |
| Amphimedon queenslandica | 13397 | 166.679 | 43615 | 95.006 | 47.162 | 3358 | 7584097 | 1612247 | V2.1 | v1 | RNAseq |
| Sycon ciliatum | 7780 | 357.509 | 32309 | 181.138 | 37.581 | 5032 | 6704019 | 2648725 | cds_v1 | v1 | NA |
| Mnemiopsis leidyi | 5100 | 155.865 | 16548 | 91.779 | 27.589 | 1324 | 2855571 | 1133446 | ML2.2 | ML2.2 | EVM |
| Hormiphora californensis | 45 | 110.691 | 14591 | 88.901 | 29.581 | 4455 | 19111845 | 4006373 | Hcv1av93 | Hcv1av93 | manual |
| Trichoplax adhaerens | 1415 | 105.612 | 12633 | 47.299 | 20.302 | 2 | 718 | 448 | JGI_AUGUSTUS | Triad1_JGI | AUGUSTUS |
| Hoilungia hongkongensis | 669 | 87.194 | 12010 | 50.820 | 20.046 | 1 | 384 | 384 | BRAKER | Hhon_v1 | BRAKER | 


