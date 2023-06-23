# steps of genome assembly #

### assembly ###
We performed a "multiple assemblers with multiple datasets approach", similar to [Guiglielmoni et al. 2021](https://doi.org/10.1186/s12859-021-04118-3). 

A total of seven different assemblers were included in the pipeline: 

* [FLYE](https://github.com/fenderglass/Flye) by [Kolmogorov et al 2019](https://doi.org/10.1038/s41587-019-0072-8)
* [CANU](https://github.com/marbl/canu) by [Koren et al 2017](https://doi.org/10.1101/gr.215087.116)
* [WTDBG2](https://github.com/ruanjue/wtdbg2) by [Ruan and Li 2020](https://doi.org/10.1038/s41592-019-0669-3)
* [MINIASM](https://github.com/lh3/miniasm) by [Li 2016](https://doi.org/10.1093/bioinformatics/btw152)
* [SHASTA](https://github.com/paoloshasta/shasta) by [Shafin et al 2020](https://doi.org/10.1038/s41587-020-0503-6)
* [RA](https://github.com/lbcb-sci/ra), does not appear to be published
* [RAVEN](https://github.com/lbcb-sci/raven) by [Vaser and Sikic 2021](http://dx.doi.org/10.1038/s43588-021-00073-4)

Error correction was done using:

* [LoRDEC](http://atgc.lirmm.fr/lordec/) by [Salmela et al 2014](https://doi.org/10.1093/bioinformatics/btu538)

Additional polishing steps with:

* [Medaka](https://github.com/nanoporetech/medaka) by [Oxford Nanopore Technologies](https://nanoporetech.com/)
* [HyPo](https://github.com/kensung-lab/hypo) Hybrid Polisher by [Kundu et al 2019 preprint](https://www.biorxiv.org/content/10.1101/2019.12.19.882506v1)

### scaffolding ###
An iterative scaffolding was performed:

* all paired-end RNAseq reads using [P_RNA_SCAFFOLDER](https://github.com/CAFS-bioinformatics/P_RNA_scaffolder) by [Zhu et al 2018](https://doi.org/10.1186/s12864-018-4567-3)
* the TransPI reference transcriptome using [L_RNA_SCAFFOLDER](https://github.com/CAFS-bioinformatics/L_RNA_scaffolder) by [Xue et al 2013](https://doi.org/10.1186/1471-2164-14-604)
* corrected long cDNA reads (full set) using `L_RNA_SCAFFOLDER`

### HiC scaffolding ###
Hi-C reads were mapped to the draft assembly using:

* [bowtie2 v2.3.5.1](https://github.com/BenLangmead/bowtie2) by [Langmead and Salzberg 2012](https://doi.org/10.1038/nmeth.1923)
* [hicstuff v2.3.0](https://github.com/koszullab/hicstuff) by [Matthey-Doret et al 2020](https://doi.org/10.1038/s41467-020-19562-7) with parameters `-e DpnII,HinfI --iterative`.

The assembly was then scaffolded using [instaGRAAL v0.1.6](https://github.com/koszullab/instaGRAAL) no-opengl branch by [Baudry et al 2020](https://doi.org/10.1186/s13059-020-02041-z) with parameters `-l 4 -n 50 -c 1 -N 5` and automatically curated using `instaGRAAL-polish`.

### assembly evaluation ###
Using code from [WRF's lavaLampPlot repo](https://github.com/wrf/lavaLampPlot) to get GC% and coverage, for plotting.

```
~/git/lavaLampPlot/hits_to_coverage.py -f v1.29/Avas.v1.29_scaffolds.fasta -g Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.PE-DNA_hisat2.bg -l 1 > Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.PE-DNA_hisat2.gc_cov.tab
# 80320060 total bp, average coverage is 101.6
Rscript ~/git/lavaLampPlot/contig_gc_coverage.R Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.PE-DNA_hisat2.gc_cov.tab 400

~/git/lavaLampPlot/hits_to_coverage.py -f v1.29/Avas.v1.29_scaffolds.fasta -g Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.PE-DNA_bowtie2.sorted.bg -l 1 > Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.PE-DNA_bowtie2.gc_cov.tab
# 80320060 total bp, average coverage is 99.7
Rscript ~/git/lavaLampPlot/contig_gc_coverage.R Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.PE-DNA_bowtie2.gc_cov.tab

~/git/lavaLampPlot/hits_to_coverage.py -f v1.29/Avas.v1.29_scaffolds.fasta -g Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.ONT-DNA_minimap2.sorted.bg -l 1 > Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.ONT-DNA_minimap2.gc_cov.tab
# 80320060 total bp, average coverage is 105.9
Rscript ~/git/lavaLampPlot/contig_gc_coverage.R Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_softmasked.ONT-DNA_minimap2.gc_cov.tab
```

