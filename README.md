# [The genome of the reef-building glass sponge *Aphrocallistes vastus* provides insights into silica biomineralization](https://royalsocietypublishing.org/doi/10.1098/rsos.230423) #
Well-annotated and contiguous genomes are an indispensable resource for understanding the evolution, development, and metabolic capacities of organisms. Sponges, an ecologically important non-bilaterian group of primarily filter-feeding sessile aquatic organisms, are underrepresented with respect to available genome resources. Here we provide a high-quality and well-annotated genome of *Aphrocallistes vastus*, a glass sponge (Porifera: Hexactinellida) that forms large reef structures off the coast of British Columbia (Canada). We show that its genome is ~80 MB, small compared to most other metazoans, and contains nearly 2,500 nested genes, more than other genomes. Hexactinellida is characterized by a unique skeletal architecture made of amorphous silicon dioxide (SiO2), and we identified 419 differentially expressed genes between the osculum, i.e., the vertical growth zone of the sponge, and the main body. Among the upregulated ones, mineralization-related genes such as glassin, as well as collagens and actins, dominate the expression profile during growth. Silicateins, suggested being involved in silica mineralization, especially in demosponges, were not found at all in the A. vastus genome and suggests that the underlying mechanisms of SiO2 deposition in the Silicea *sensu stricto* (Hexactinellida+Demospongiae) may not be homologous.

## Description ##
This repository includes all the code used to analyze the genome of *Aphrocallistes vastus* for the article [The genome of the reef-building glass sponge *Aphrocallistes vastus* provides insights into silica biomineralization](https://royalsocietypublishing.org/doi/10.1098/rsos.230423). All data, figures, and supplemental files can be found in this repo, as well as in the [Zenodo link](https://doi.org/10.5281/zenodo.7970685). Raw reads for assembly and annotation can be found at [ENA Project PRJEB61987](https://www.ebi.ac.uk/ena/browser/view/PRJEB61987).

**Most relevant files for download are in the [Releases tab](https://github.com/PalMuc/Aphrocallistes_vastus_genome/releases), including the assembly, gene annotation as GFF, transcripts, and proteins.** These files can also be found individually in the [assembly folder](https://github.com/PalMuc/Aphrocallistes_vastus_genome/tree/main/assembly).

**Predicted proteins of the 14 additional Hexactinellid transcriptomes are in the [transciptomes_SONNE folder](https://github.com/PalMuc/Aphrocallistes_vastus_genome/tree/main/transciptomes_SONNE).**

### Other folders include: ###

* [genomic_overview](https://github.com/PalMuc/Aphrocallistes_vastus_genome/tree/main/genomic_overview), data for Table 1, Figure 1, part of Figure 2, SFigs 1-3
* [synteny](https://github.com/PalMuc/Aphrocallistes_vastus_genome/tree/main/synteny), data related to the synteny analysis in Figure 2 and SFig 5
* [biomineralization](https://github.com/PalMuc/Aphrocallistes_vastus_genome/tree/main/biomineralization), data, code, and gene trees, used to make parts of Figures 3, 4, and 5
* [differential_gene_expression](https://github.com/PalMuc/Aphrocallistes_vastus_genome/tree/main/differential_gene_expression), data and code for DEG analysis, used in Figures 3-6, and SFigures 8-9
* [ortholog_clusters](https://github.com/PalMuc/Aphrocallistes_vastus_genome/tree/main/ortholog_clusters), data for generating gene clusters
* [hexa_phylogeny](https://github.com/PalMuc/Aphrocallistes_vastus_genome/tree/main/hexa_phylogeny), data for Hexactinellid phylogeny in Figure 1
* [nested_genes](https://github.com/PalMuc/Aphrocallistes_vastus_genome/tree/main/nested_genes), data for nested gene analysis in Table 1
* [trans_splicing_leader](https://github.com/PalMuc/Aphrocallistes_vastus_genome/tree/main/trans_splicing_leader), data for trans-splicing in SFig 4
  
* [figures_for_paper](https://github.com/PalMuc/Aphrocallistes_vastus_genome/tree/main/figures_for_paper), PDF versions of all figures and subparts
* [supplements_for_paper](https://github.com/PalMuc/Aphrocallistes_vastus_genome/tree/main/supplements_for_paper), all supplemental figures (SFig 01-SFig 13), tables (STab 01 - STab 10), two trees, and three alignments


### Full citation ###
Francis, W.R.; Eitel, M.; Vargas, S.; Garcia-Escudero, C. A.; Conci, N.; Deister, F.; Mah, J. L.; Guiglielmoni, N., Krebs, S.; Blum, H.; Leys, S. P.; Wörheide, G. (2023) The genome of the reef-building glass sponge *Aphrocallistes vastus* provides insights into silica biomineralization. Royal Society Open Science, 10(6), 230423. DOI: [10.1098/rsos.230423](https://royalsocietypublishing.org/doi/10.1098/rsos.230423)

Francis, Warren R. 1,§; Eitel, Michael 1,§; Vargas, Sergio 1; Garcia-Escudero, Catalina, A. 1; Conci, Nicola 1; Deister, Fabian 1; Mah, Jasmine L. 2; Guiglielmoni, Nadège 3; Krebs, Stefan 4; Blum, Helmut 4; Leys, Sally P. 2; Wörheide, Gert 1,5,6,%

1. Department of Earth and Environmental Sciences, Paleontology and Geobiology, Ludwig-Maximilians-Universität München, Munich, Germany.
2. Department of Biological Sciences, University of Alberta, Edmonton, AB T6G 2E9, Canada.
3. Service Evolution Biologique et Ecologie, Université libre de Bruxelles (ULB), 1050 Brussels, Belgium
4. Laboratory for Functional Genome Analysis (LAFUGA), Gene Center, Ludwig-Maximilians-Universität München, Munich, Germany
5. GeoBio-Center, Ludwig-Maximilians-Universität München, Munich, Germany
6. Staatliche Naturwissenschaftliche Sammlungen Bayerns (SNSB)–Bayerische Staatssammlung für Paläontologie und Geologie, Munich, Germany

§ joint first authors

% Corresponding author

This work was supported by the European Union's Horizon 2020 research and innovation programme under Marie Skłodowska-Curie grant agreement no. 764840 (ITN IGNITE) and through the LMU Munich's Institutional Strategy LMUexcellent within the framework of the German Excellence Initiative to G.W., as well as through a NSERC Discovery grant (grant no. 2016-05446) to S.P.L.

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.
