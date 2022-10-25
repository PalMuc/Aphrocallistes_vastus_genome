# synteny analysis of *A. vastus* vs others #

Using [diamond](https://github.com/bbuchfink/diamond) and [WRF's genomeGTFtools](https://github.com/wrf/genomeGTFtools):

### compared with *Oopsacas minuta* ###
This is the NCBI-uploaded version of the *O. minuta* [JAKMXF01](https://www.ncbi.nlm.nih.gov/Traces/wgs/JAKMXF01), which should be the [preprint version by Santini et al](https://www.biorxiv.org/content/10.1101/2022.07.26.501511v1.full). Some additional processing steps from the GenBank format were done, following [WRF's annotation instructions](https://bitbucket.org/wrf/genome-reannotations/src/master/jbrowse-tracks/oopsacas/).



### with *Ephydatia muelleri* ###
This is the [version v1](https://spaces.facsci.ualberta.ca/ephybase/) from the *Ephydatia muelleri* genome paper by [Kenny et al 2020](https://doi.org/10.1038/s41467-020-17397-w)

```
~/diamond-v2.0.13/diamond-v2.0.13 blastp -q v1.29/Avas.v1.29_annotations.prot.fr.fasta -d ~/genomes/ephydatia_muelleri/augustus_sysnames_prots.fasta.dmnd -o Avas.v1.29_annotations.prot.blastp_vs_Emue1.tab

~/git/genomeGTFtools/scaffold_synteny.py -b Avas.v1.29_annotations.prot.blastp_vs_Emue1.tab -q v1.29/Avas.v1.29_annotations.fr.gff -d ~/genomes/ephydatia_muelleri/augustus/Emu_augustus_vs_nb_sysnames_gene_only.gff -f v1.29/Avas.v1.29_scaffolds.fasta -F ~/genomes/ephydatia_muelleri/ASM_HIC_394/renumbered_final_nb.fa -l 80 -L 250 -G 50 > avas1.29_vs_Emue1_scaffold2d_points.tab

Rscript ~/git/genomeGTFtools/synteny_2d_plot.R avas1.29_vs_Emue1_scaffold2d_points.tab Aphrocallistes Ephydatia 140

~/git/genomeGTFtools/scaffold_synteny.py -b Avas.v1.29_annotations.prot.blastp_vs_Emue1.tab -q v1.29/Avas.v1.29_annotations.fr.gff -d ~/genomes/ephydatia_muelleri/augustus/Emu_augustus_vs_nb_sysnames_gene_only.gff -f v1.29/Avas.v1.29_scaffolds.fasta -F ~/genomes/ephydatia_muelleri/ASM_HIC_394/renumbered_final_nb.fa -l 80 -L 250 -G 50 -R > avas1.29_vs_Emue1_scaffold2d_points.rand.tab

Rscript ~/git/genomeGTFtools/synteny_2d_plot.R avas1.29_vs_Emue1_scaffold2d_points.rand.tab Aphrocallistes Ephydatia 160

Rscript ~/git/genomeGTFtools/synteny_2d_fishers_test.R avas1.29_vs_Emue1_scaffold2d_points.tab A.vastus E.muelleri

```

