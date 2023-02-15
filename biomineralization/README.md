# steps for analyzing biomineralization genes #
Mineralizing proteins are often repetitive, such as the [pearl oyster aspein](https://www.uniprot.org/uniprotkb/Q76K52/entry) (from [Tsukamoto 2004](https://doi.org/10.1016/j.bbrc.2004.06.072)):

```
>Pinctada_fucata|Q76K52_PINFU Shell matrix protein OS=Pinctada fucata OX=50426 GN=aspein PE=2 SV=1
MKGIAILMCLAALVAVSVTFPVADQTTNELGSSGAAAAGAVVSEPSDAGDAADAGDADAA
DADAADADADADADADNDDGDDDDDDDDDDSGDDDSGDDDDSGDDDDSGDDDDSGDDDDS
GDDDDSGDDDGDDDSEDDDDSGDDDSGDGDDGDSGDDDDDDDSGDDDDDDDSGDDDDDSG
DDDDGDSGDDDSGDDDGDDDDSGDDDSGDDDDSGDDDSGDDDDSGDDDSGDDESGDDDSG
DDDDSGDDDSGDDDSGDDGDDDDSGDDDDSGDDDDSGDDDDDDDDSGDDDDGDSGDDDSG
DDDGDDDDSGDDDSGDDESGDDDSEDDDSGDDDSGDDDSGDDDSDSGDDDSGDDDSGDDD
SGDDDGDDGDDDADSGDDDDDDDDDDGDDGDDDSGDDDGDDSDDDDDDDDDDQ
```

### glassin ###
Plots were made using code from [WRF's genomeGTFtools repo](https://github.com/wrf/genomeGTFtools). The partial [glassin protein](https://www.ncbi.nlm.nih.gov/protein/BAS21353.1) from *Euplectella curvistellata* was identified by [Shimizu 2015](http://www.pnas.org/lookup/doi/10.1073/pnas.1506968112), and the loci shown are the complete gene in both species.

```
~/git/genomeGTFtools/extract_coordinates.py -g JAKMXF01.1.gbff.gff -b 215000 -e 240000 -s JAKMXF010000288.1 > JAKMXF01.1.gbff.scaf_288.215k_240k_glassin.tab
Rscript ~/git/genomeGTFtools/draw_annotation_blocks.R JAKMXF01.1.gbff.scaf_288.215k_240k_glassin.tab

~/git/genomeGTFtools/extract_coordinates.py -g Avas.v1.29_annotations.fr.gff -b 2465000 -e 2490000 -s Aphrocallistes_vastus_HiC-scaffold_014 > Avas.s014.2465k_2490k_glassin.tab
Rscript ~/git/genomeGTFtools/draw_annotation_blocks.R Avas.s014.2465k_2490k_glassin.tab
```

Molecular weight of sequences and potential cleavage products assuming a cut at `XG'G`, where `X` is non-polar, and `'` is the cutting site:

| sequence ID                   | molecular weight | pI    |
| :---                                 |     ---: |   ---: |
| Euplectella_curvistellata_LC012028.1 | 14772.07 | 10.140 |
| Euplectella_curvistellata_BAS21356.1 | 19437.71 | 6.480 |
| Euplectella_curvistellata_LC010923.1 | 26369.67 | 6.160 |
| Euplectella_curvistellata_BAS21353.1 | 32797.34 | 4.890 |
| Euplectella_curvistellata_LC012024.1 | 32797.34 | 4.890 |
| Oopsacas_minuta_LOD99_3750.mRNA.1 | 47800.08 | 5.530 |
| O.minuta_glassin-cut | 25692.17 | 5.890 |
| Aphrocallistes_vastus_Avas.s014.g618.i1 | 48056.01 | 5.700 |
| A.vastus_glassin-cut | 25754.69 | 6.390 |
| Rosella_fibulata_TR13145_c0_g2_i1_2 | 50869.66 | 6.110 |
| R.fibulata_glassin-cut | 28446.33 | 6.340 |

Gene appears to be single copy in Avas. Only hit in the genome is the known locus on scaffold 014:

```
tblastn -query glassin_euplectella_only.fasta -db ../../v1.29/Avas.v1.29_scaffolds.fasta -outfmt 6 -evalue 1e-1
Euplectella_curvistellata|BAS21353.1	Aphrocallistes_vastus_HiC-scaffold_014	40.116	172	102	1	121	292	2478578	2478066	2.18e-36	139
Euplectella_curvistellata|LC012024.1	Aphrocallistes_vastus_HiC-scaffold_014	40.116	172	102	1	121	292	2478578	2478066	2.18e-36	139
```

### TIL domain protein ###
The identification of proteins with compositional bias was done with the script [protein_repeat_finder.py](https://bitbucket.org/wrf/sequences/src/master/protein_repeat_finder.py), making the table output only. Option `-m` restricts proteins with at least 50AAs, avoiding compositional bias flags for very short proteins.

`protein_repeat_finder.py Avas.v1.29_annotations.prot.fr.fasta -m 50 > Avas.v1.29_annotations.prot.fr.rep_flags_m50.tab`

Then include the protein sequences with the flags:

`protein_repeat_finder.py -s Avas.v1.29_annotations.prot.fr.fasta -m 50 > Avas.v1.29_annotations.prot.fr.rep_flags_m50.fasta`

```
>Avas.s003.g337.i1 0 3698 Com:H-15.09
MKIFGEIVLLLSILLLANAARIPLCSNDEIYSECERCPLMCENLLSDDICGTECMKGCFC
KEGLARNRTGYCVDVFTCPLEYKQECEGENEEFKECGSRCSTTCENGVEIEICTKGCYPG
CMCKEGFVRNSAGVCAALAECGVSGVSQPDTECGENQVQTECMTCPETCDNRETKSDCLL
ECRSGCACAPGYILASNITNECINISKCQYKHKHKHGHRHNHQHERNHHSHQHGRNHSHQ
HERNHEHQHERNHEHQHEHNHEHQHERNHEHQHEHNHEHQHEHNHEHQHEHGQDSHEHQH
ENQHGNQHENQHENHHGNQHEGHEHGSSHGSHPRRHNSHSHTQSQGNSNSPDDSHSNKES
NHRRHDKHDRS
```

Schematic of the loci is made by:

```
~/git/genomeGTFtools/extract_coordinates.py -g JAKMXF01.1.gbff.gff -b 410000 -e 430000 -s JAKMXF010000144.1 > JAKMXF01.1.gbff.scaf_144.410k_430k_TIL.tab
Rscript ~/git/genomeGTFtools/draw_annotation_blocks.R JAKMXF01.1.gbff.scaf_144.410k_430k_TIL.tab

~/git/genomeGTFtools/extract_coordinates.py -g Avas.v1.29_annotations.fr.gff -b 875000 -e 895000 -s Aphrocallistes_vastus_HiC-scaffold_003 > Avas.s003.875k_895k_TIL.tab
Rscript ~/git/genomeGTFtools/draw_annotation_blocks.R Avas.s003.875k_895k_TIL.tab
```

In other proteins, the disulfide bonds found were: C5–C38, C15–C33, C18–C29, C22–C60 and C40–C54.

### carbonic anhydrase ###
For schematic of loci:

```
~/git/genomeGTFtools/extract_coordinates.py -g JAKMXF01.1.gbff.gff -b 75000 -e 125000 -s JAKMXF010000310.1 -G | sed s/".mRNA.1"//g | sed s/"LOD99_"//g > JAKMXF01.1.gbff.scaf_310.75k_125k_CA.tab
Rscript ~/git/genomeGTFtools/draw_annotation_blocks.R JAKMXF01.1.gbff.scaf_310.75k_125k_CA.tab
~/git/genomeGTFtools/extract_coordinates.py -g JAKMXF01.1.gbff.gff -b 1375000 -e 1425000 -s JAKMXF010000222.1 -G | sed s/".mRNA.1"//g | sed s/"LOD99_"//g > JAKMXF01.1.gbff.scaf_222.1375k_1425k_CA.tab
Rscript ~/git/genomeGTFtools/draw_annotation_blocks.R JAKMXF01.1.gbff.scaf_222.1375k_1425k_CA.tab

~/git/genomeGTFtools/extract_coordinates.py -g Avas.v1.29_annotations.fr.gff -b 5325000 -e 5375000 -s Aphrocallistes_vastus_HiC-scaffold_001 -G | sed s/"Avas.s001."//g | sed s/".i1"//g > Avas.s001.5325k_5375k_CA.tab
Rscript ~/git/genomeGTFtools/draw_annotation_blocks.R Avas.s001.5325k_5375k_CA.tab
```

For supplemental figure alignment with [MAFFT](https://mafft.cbrc.jp/alignment/software/) and tree with [iqtree2](http://www.iqtree.org/):

```
mafft-linsi silicases_and_cah.fasta > silicases_and_cah.aln
~/iqtree-2.2.0-Linux/bin/iqtree -s silicases_and_cah.aln -mset WAG,LG -b 100 -nt 4
```

The structure with sequence identity of carbonic anhydrase-silicase is made using Alphafold structure for [Aqu2.1.41692_001](https://alphafold.ebi.ac.uk/entry/A6QR76), which matches the UniProt protein [A6QR76](https://www.uniprot.org/uniprotkb/A6QR76/entry).

`~/git/pdbcolor/pdb_site_identity.py -a demo_subtree.aln -p AF-A6QR76-F1-model_v4.pdb -s AQUE_Aqu2.1.41692_001 -w AF-A6QR76-F1-model_v4.color_by_identity.pml --base-color blue --force-recode`

This generates a [PyMol script](https://pymol.org/2/), run as:

`@/home/dummy/git/Aphrocallistes_vastus_genome/biomineralization/AF-A6QR76-F1-model_v4.color_by_identity.pml`


### actins ###
Matching the actin sequence from [Ehrlich 2022](https://onlinelibrary.wiley.com/doi/full/10.1002/advs.202105059), three nearly identical actins are found in Avas.

```
# successive identification of exact motifs
findmotif.py -m AGFAGDDAPRAVFPSIVGRPR paralogs_holozoa_allprots_v10_00016.fasta > round1
findmotif.py -m VAPEEHPVLLTEAPLNPK round1 > round2
findmotif.py -m GYSFTTTAER round2 > round3
getAinB.py Avas round3 -s > actins_matching_ehrlich2022.fasta
```


### collagen loci ###
Schematic is generated using:

```
~/git/genomeGTFtools/extract_coordinates.py -g Avas.v1.29_annotations.fr.gff -b 3349000 -e 3400000 -s Aphrocallistes_vastus_HiC-scaffold_001 -G --make-linear | sed s/"Avas.s001."//g | sed s/".i1"//g > Avas.s001.3349k_3400k_col.tab
Rscript ~/git/genomeGTFtools/draw_annotation_blocks.R Avas.s001.3349k_3400k_col.tab
~/git/genomeGTFtools/extract_coordinates.py -g Avas.v1.29_annotations.fr.gff -b 3860000 -e 3910000 -s Aphrocallistes_vastus_HiC-scaffold_006 -G --make-linear | sed s/"Avas.s006."//g | sed s/".i1"//g > Avas.s006.3860k_3910k_col.tab
Rscript ~/git/genomeGTFtools/draw_annotation_blocks.R Avas.s006.3860k_3910k_col.tab
```

### DUOX loci ###
Back-to-back loci in several genomes, schematic is generated using:

```
~/git/genomeGTFtools/extract_coordinates.py -g Avas.v1.29_annotations.fr.gff -b 968000 -e 1008000 -s Aphrocallistes_vastus_HiC-scaffold_019 -G --make-linear | sed s/"Avas.s019."//g | sed s/".i1"//g > duox/Avas.s019.968k_1008k_duox.tab
Rscript ~/git/genomeGTFtools/draw_annotation_blocks.R duox/Avas.s019.968k_1008k_duox.tab
~/git/genomeGTFtools/extract_coordinates.py -g JAKMXF01.1.gbff.gff -b 111000 -e 151000 -s JAKMXF010000322.1 -G --make-linear | sed s/".mRNA.1"//g | sed s/"LOD99_"//g > duox/JAKMXF01.1.gbff.scaf_322.111k_151k_duox.tab
Rscript ~/git/genomeGTFtools/draw_annotation_blocks.R duox/JAKMXF01.1.gbff.scaf_322.111k_151k_duox.tab

~/git/genomeGTFtools/extract_coordinates.py -g ~/genomes/ephydatia_muelleri/Emu_augustus_vs_nb_sysnames.gff -b 2664000 -e 2704000 -s scaffold_0023 -G --make-linear | sed s/".t1"//g > duox/Emu_augustus_vs_nb_sysnames.scaffold_0023.2664k_2704k_duox.tab
Rscript ~/git/genomeGTFtools/draw_annotation_blocks.R duox/Emu_augustus_vs_nb_sysnames.scaffold_0023.2664k_2704k_duox.tab
~/git/genomeGTFtools/extract_coordinates.py -g ~/genomes/amphimedon_queenslandica_PORI/Amphimedon_queenslandica.Aqu1.45.cleaned.gff -b 1 -e 40000 -s Contig13501 -G --make-linear | sed s/"transcript:Aqu2.1."//g > duox/Amphimedon_queenslandica.Aqu1.45.Contig13501.0_40k_duox.tab
Rscript ~/git/genomeGTFtools/draw_annotation_blocks.R duox/Amphimedon_queenslandica.Aqu1.45.Contig13501.0_40k_duox.tab

~/git/genomeGTFtools/extract_coordinates.py -g ~/genomes/salpingoeca_rosetta_CHOA/Salpingoeca_rosetta.Proterospongia_sp_ATCC50818.36.gff3 -b 1406000 -e 1446000 -s supercont1.15 -G --make-linear | sed s/"transcript:"//g > duox/Salpingoeca_rosetta.supercont1.15.1406k_1446k_duox.tab
Rscript ~/git/genomeGTFtools/draw_annotation_blocks.R duox/Salpingoeca_rosetta.supercont1.15.1406k_1446k_duox.tab
```

Check [Sebe-Pedros 2018](http://www.nature.com/articles/s41559-018-0575-6) data for *Amphimedon*, finding clusters 23, 41-46, 52 in adult containing both genes, and clusters 15-17 in larvae. These are adult pinacocytes according to Figure 2G of that paper, and "non-ciliated epithelium" cells in larvae.

```
grep Aqu2.1.39906_001 *.csv
41559_2018_575_MOESM5_ESM_aque_adult_clusters.C23.csv:Aqu2.1.39906_001;8400;593;7.0595238095;2.9800952278;An_peroxidase/An_peroxidase/EF-hand_7/Ferric_reduct/FAD_binding_8/NAD_binding_6/; cDNA FLJ76088, highly similar to Homo sapiens dual oxidase 1 (DUOX1), transcript variant 1, mRNA; Dual oxidase
41559_2018_575_MOESM5_ESM_aque_adult_clusters.C41.csv:Aqu2.1.39906_001;8400;1300;15.4761904762;2.9901912183;An_peroxidase/An_peroxidase/EF-hand_7/Ferric_reduct/FAD_binding_8/NAD_binding_6/; cDNA FLJ76088, highly similar to Homo sapiens dual oxidase 1 (DUOX1), transcript variant 1, mRNA; Dual oxidase
41559_2018_575_MOESM5_ESM_aque_adult_clusters.C42.csv:Aqu2.1.39906_001;8400;758;9.0238095238;3.9833773117;An_peroxidase/An_peroxidase/EF-hand_7/Ferric_reduct/FAD_binding_8/NAD_binding_6/; cDNA FLJ76088, highly similar to Homo sapiens dual oxidase 1 (DUOX1), transcript variant 1, mRNA; Dual oxidase
41559_2018_575_MOESM5_ESM_aque_adult_clusters.C43.csv:Aqu2.1.39906_001;8400;214;2.5476190476;1.586187134;An_peroxidase/An_peroxidase/EF-hand_7/Ferric_reduct/FAD_binding_8/NAD_binding_6/; cDNA FLJ76088, highly similar to Homo sapiens dual oxidase 1 (DUOX1), transcript variant 1, mRNA; Dual oxidase
41559_2018_575_MOESM5_ESM_aque_adult_clusters.C44.csv:Aqu2.1.39906_001;8400;806;9.5952380952;3.0915097109;An_peroxidase/An_peroxidase/EF-hand_7/Ferric_reduct/FAD_binding_8/NAD_binding_6/; cDNA FLJ76088, highly similar to Homo sapiens dual oxidase 1 (DUOX1), transcript variant 1, mRNA; Dual oxidase
41559_2018_575_MOESM5_ESM_aque_adult_clusters.C45.csv:Aqu2.1.39906_001;8400;500;5.9523809524;2.6965399559;An_peroxidase/An_peroxidase/EF-hand_7/Ferric_reduct/FAD_binding_8/NAD_binding_6/; cDNA FLJ76088, highly similar to Homo sapiens dual oxidase 1 (DUOX1), transcript variant 1, mRNA; Dual oxidase
41559_2018_575_MOESM5_ESM_aque_adult_clusters.C46.csv:Aqu2.1.39906_001;8400;437;5.2023809524;2.3642103585;An_peroxidase/An_peroxidase/EF-hand_7/Ferric_reduct/FAD_binding_8/NAD_binding_6/; cDNA FLJ76088, highly similar to Homo sapiens dual oxidase 1 (DUOX1), transcript variant 1, mRNA; Dual oxidase
41559_2018_575_MOESM5_ESM_aque_adult_clusters.C51.csv:Aqu2.1.39906_001;8400;574;6.8333333333;1.49058811;An_peroxidase/An_peroxidase/EF-hand_7/Ferric_reduct/FAD_binding_8/NAD_binding_6/; cDNA FLJ76088, highly similar to Homo sapiens dual oxidase 1 (DUOX1), transcript variant 1, mRNA; Dual oxidase
41559_2018_575_MOESM5_ESM_aque_adult_clusters.C52.csv:Aqu2.1.39906_001;8400;277;3.2976190476;1.4118036705;An_peroxidase/An_peroxidase/EF-hand_7/Ferric_reduct/FAD_binding_8/NAD_binding_6/; cDNA FLJ76088, highly similar to Homo sapiens dual oxidase 1 (DUOX1), transcript variant 1, mRNA; Dual oxidase

grep Aqu2.1.39905_001 *.csv
41559_2018_575_MOESM5_ESM_aque_adult_clusters.C22.csv:Aqu2.1.39905_001;4350;80;1.8390804598;1.1792683627;DuoxA/; Dual oxidase maturation factor 1 alpha; LD27791p
41559_2018_575_MOESM5_ESM_aque_adult_clusters.C23.csv:Aqu2.1.39905_001;4350;302;6.9425287356;2.7798331381;DuoxA/; Dual oxidase maturation factor 1 alpha; LD27791p
41559_2018_575_MOESM5_ESM_aque_adult_clusters.C41.csv:Aqu2.1.39905_001;4350;690;15.8620689655;2.809069348;DuoxA/; Dual oxidase maturation factor 1 alpha; LD27791p
41559_2018_575_MOESM5_ESM_aque_adult_clusters.C42.csv:Aqu2.1.39905_001;4350;379;8.7126436782;3.6715705566;DuoxA/; Dual oxidase maturation factor 1 alpha; LD27791p
41559_2018_575_MOESM5_ESM_aque_adult_clusters.C43.csv:Aqu2.1.39905_001;4350;128;2.9425287356;1.5341213953;DuoxA/; Dual oxidase maturation factor 1 alpha; LD27791p
41559_2018_575_MOESM5_ESM_aque_adult_clusters.C44.csv:Aqu2.1.39905_001;4350;341;7.8390804598;2.4705742039;DuoxA/; Dual oxidase maturation factor 1 alpha; LD27791p
41559_2018_575_MOESM5_ESM_aque_adult_clusters.C45.csv:Aqu2.1.39905_001;4350;271;6.2298850575;2.4221640225;DuoxA/; Dual oxidase maturation factor 1 alpha; LD27791p
41559_2018_575_MOESM5_ESM_aque_adult_clusters.C46.csv:Aqu2.1.39905_001;4350;220;5.0574712644;2.0648894919;DuoxA/; Dual oxidase maturation factor 1 alpha; LD27791p
41559_2018_575_MOESM5_ESM_aque_adult_clusters.C52.csv:Aqu2.1.39905_001;4350;155;3.5632183908;1.3097754563;DuoxA/; Dual oxidase maturation factor 1 alpha; LD27791p
```

The DUOX gene in *Spongilla lacustris* was specifically highlighted in a main figure 1E in [Musser 2021](https://zenodo.org/record/5094890), with high expression in "myopeptidocytes".

### silicateins in Spongilla and other genes ###
Dataset [on Zenodo](https://zenodo.org/record/5094890), file `spongilla_transcriptome_annotated_200501.fasta` contains 62180 transcripts, though only 24055 could map to the [Ephydatia genome](https://spaces.facsci.ualberta.ca/ephybase/)

```
~/minimap2-2.23_x64-linux/minimap2 -a -x splice --secondary=no renumbered_final_nb.fa spongilla_transcriptome_annotated_200501.fasta | ~/samtools-1.14/samtools sort - -o spongilla_transcriptome_annotated_200501.vs_emue1.bam
~/samtools-1.14/samtools flagstat spongilla_transcriptome_annotated_200501.vs_emue1.bam
#62739 + 0 in total (QC-passed reads + QC-failed reads)
#62180 + 0 primary
#559 + 0 supplementary
#24055 + 0 mapped (38.34% : N/A)
#23496 + 0 primary mapped (37.79% : N/A)
~/git/pinfish/spliced_bam2gff/spliced_bam2gff -M -s spongilla_transcriptome_annotated_200501.vs_emue1.bam > spongilla_transcriptome_annotated_200501.vs_emue1.gtf
~/git/genomeGTFtools/misc/stringtie_gtf_to_gff3.py spongilla_transcriptome_annotated_200501.vs_emue1.gtf > spongilla_transcriptome_annotated_200501.vs_emue1.gff
```

To examine sclerocyte specific genes, and proteins for repeats or compositional bias. None of the below candidates are especially promising. Potentially the gene has low expression or was incomplete, and the genome would be needed.

```
grep "Scl" Data_S1_sheet3_diff_exp_cell_types.tab > Data_S1_sheet3_diff_exp_cell_types.scl_only.tab
cut -f 2 Data_S1_sheet3_diff_exp_cell_types.scl_only.tab | cut -f 1 -d " " > Data_S1_sheet3_diff_exp_cell_types.scl_only.names
getAinB.py Data_S1_sheet3_diff_exp_cell_types.scl_only.names ../master/transcriptome_proteome_and_phylome/spongilla_transcriptome_annotated_200501.fasta.gz -s > Data_S1_sheet3_diff_exp_cell_types.scl_only.nucl.fasta
prottrans.py -m Data_S1_sheet3_diff_exp_cell_types.scl_only.nucl.fasta > Data_S1_sheet3_diff_exp_cell_types.scl_only.prot.fasta
protein_repeat_finder.py -m 50 -s Data_S1_sheet3_diff_exp_cell_types.scl_only.prot.fasta > Data_S1_sheet3_diff_exp_cell_types.scl_only.prot.fasta.repeats
```

The transcripts in both species are missing start codons and stop codons, so not likely to be read in the sense with the repeats.

```
>c93746_g1_i1|c93746_g1 1 551 Com:S-21.09
TLQCTCTECYHNHLLWTTCTEHEHNQYSVKISTHYSDHSSSSSSSPPAASSSSSSLCCIT
LERMARALASAPAIPSRSSGCAQDSRESSVETVFSFIPSLMLDSKLLRASGTAGVGGGAC
LREGGGGA
>Em_i100k_d5.40034.1 5
FVLDEFINDSLNPARDDTGACSLKATHKNTISWLLTVKVSHTHYSDHSSSSSSSPPASSS
SSSLCCMTLLRMARALASAPAIPSRSSEGAHDSRESSTDTVFSCIPALMLDSKLLRESGT
GGAAADAWCLPGGGG
```

Expression level high outside of sclerocytes, identified as Zinc-finger-domain-containing

```
>c98317_g1_i1|c98317_g1 4 1037 Com:S-28.22
MASVSHTTVNSETARRYLKGKQTQQDASSYRCQKCLQYGHYTYECTGKRKYLYRPSRTKN
MKNNKSETPEVSRSTISKKERSKEEPKKKRARRASTTSESDESSSSSSDTSSSDSDSSDS
ESDSTASDSSSESSSSSESDSSTSSSSDSH
>Em0022g158a
MSSVTHTSVSSETARRYVKGKQTQQDARSYRCQKCLQYGHYTYECTGKRKYLYRPSRTKN
MKKGNTDTSEVSRSANASKKERSREEPKKKRARRASTSSGSDSSSSSSDTSSSDSESSES
ESDSDTTSDSSSSSSSESQDSTESDSSTSSSSDSD
```

Identified as [calreticulin](https://www.uniprot.org/uniprotkb/P27797/entry), very high expression in *E. muelleri*, high expression outside of sclerocytes

```
>c100212_g2_i1|c100212_g2 5 1363 Rep:D-21 Com:D-17.41
MERIVVFLALIVVALAEPTIYFKETFEDATWKDRWVQSKHKSDFGEFVWSAGKFYGDA
ELDKGIKTSQDARFYALSAGFDKFSNEGKTLVIQFTVKHEQNIDCGGGYVKLFPSTQDQK
DLNGESPYNIMFGPDICGPGTKKVHVIFNYKEKNLLVKKEIRCKDDEYTHLYTLIVRPDN
TYEVKIDNSKVESGTLESDFDFLPAKTIPDPAASKPSDWDDKAKIDDPSDTKPAEWDKPE
YIADPDAKKPEDWDDDTDGEWEPPNVPNPAYKGEWKPKQIDNPNYKGVWVHPQVPNPAYQ
PDDGLYSYDDFGVVGFDLWQVKSGTIFDNVLITDDVEYAETFGADTWGKTKDAEKEMKKK
LDEKQREEDEARRKAEEAAKAKDDDNDDDDDDDDDDDDDDDDDDDDDHDGHDHDHDGDAA
PKEEL
>Em0011g774a
MERIVIVLLVLLGVALADPTIYFKETFEDAAWKDRWVQSKHKSDFGEFTWSAGKFYGDAE
KDKGIKTSQDARFYALSAGFDKFSNEDKTLVIQFTVKHEQSIDCGGGYVKLFPSTVDKKD
VHGESPYNIMFGPDICGPGTKKVHVIFSYKGKNLLVKKEIRCKDDEYTHLYTLIVRPDNT
YEVKIDNSKVESGSLENDFDFLPAKTIPDPAASKPSDWDDRAKIDDPSDTKPEDWDKPEY
IPDPDAKKPEDWDDDTDGEWEPPNVPNPAYKGEWKPKQIDNPKYQGVWVHPQVPNPEYQP
DDSLYSYDDFGSIGFDLWQVKSGTIFDNVLITDDVEYAETFGADTWGKTKDAEKEMKKKL
DEEQRKEEEAKRKAEEASKAKDDDDDDDDDDDDDDNDDHDHDHEAEAAPKEEL
```

Identified as SLC39A7 solute carrier family 39 member 7 similar to zinc transporter, high expression in sclerocytes, low-medium expression outside. Histidine repeats are also present in human [S39A7_HUMAN](https://www.uniprot.org/uniprotkb/Q92504/entry)

```
>c103758_g1_i1|c103758_g1 3 2301 Rep:L-6 Rep:H-6 Com:H-10.93
MQGAKRRNYLTANMHVLLLLLLCLPAVLCHGDHSHDHGGHSHDHDHGG
HGLHSHFRSKKSHSHDHAHGHDHGHSHDHGHNHKHDHVAQKGHEHGQGGTPSSGHGTHVS
KSAPLENSWEVWGYALGSTALISAAPFVILFFIPITDATEHSTLLKVLLSFASGGLLGDA
FLHLIPHAVSPHHHHHHDGEEGDHVTSHDHDEHGGDAHNHMQDMVVGLWVLSGIIAFLLV
EKFVRLAKGSGHSHGCHAHGEPVTKKDNTAEGEGEEEDEKKEQTHDRKEKDGSLRKRGAS
AKEKADPKASPKASPKASTKAGKDTTTPAPPPPTTSDIKVAGYLNLAADCAHNFTDGLAI
GASYLAGHSIGVVTTLTILFHELPHEIGDFAILVQSGYPKRKAILLQLVTAVGALAGTLC
GLCAQSVGVAATAWILPFTAGGFIYIATVSVIPVLLENSSPGQSVLEILGMFAGIGMMVV
IAFIE
>Em0016g917a
MQQGAKKRGHTKTNMHAAPLALLLLLLCLPAVFGHGGRAHDHGGHKHGGHNHGHDHDHEG
HGLHSQFRTKHPHGHDHGHSHDHGHSHDHGHSHAHAHDHAHDHGHSHHHAHDHKHDHVDD
KGHGHGAVNSVPENGWEMWGYALGSTALISAAPFVILFFIPITDATTEHSTLLKLLLSFA
SGGLLGDAFLHLIPHAVSPHHHHDDGKEGGHVMSHDHEGDAHDHMQDMVVGLWVLSGIIA
FLLVEKFVRLVKGGHSHGHGGHCHGEPTTEKDEIPSSEKEEAKTERSHKEDGDDKDGTVR
KRKADSKEKADSKEKADSKEKAESKEKADSKEKAATRKPAPASSDIKVAGYLNLAADCAH
NFTDGLAIGASYLAGHSVGVVTTLTILFHELPHEVGDFAILVQSGYPKRKAILLQLVTAI
GALAGTLCGLFADSVGVAATAWILPFTAGGFIYIATVSVIPVLLENSSPGQSILEILAMF
AGIGMMVVIAFVE
```


