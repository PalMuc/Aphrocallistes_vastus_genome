# microbes associated with *A. vastus* #

### Table 1: Overview of microbial bins ###
| parameter      | bin 1            | bin 2             | *A. vastus*       |
| :---           | ---              | ---               | ---               |
| Genome size    | 2.22 Mb          | 1.57 Mb           | 80.32 Mb          |
| N scaffolds    | 21               | 1                 | 186               |
| N50            | 289 kb           | 1.57 Mb           | 3.37 Mb           |
| Largest scaffold | 0.605 Mb       | 1.57 Mb           | 8.19 Mb           |
| GC percent     | 48.5             | 56.2              | 37.4              |
| Genes (prodigal) | 2131            | 1589              | 19578            |
| w/ annotation (KEGG) | 1207 (56%)  | 987 (62%)         | -             |
| Putative group | b-proteobacteria | a-proteobacteria  | sponge            |
| Putative ID (ATPase)  | Nitrosomonas (75%) | Rhodospirillaceae (73%) | -  |
| Hit accession (ATPase) | [WP_011110849.1](https://www.ncbi.nlm.nih.gov/protein/WP_011110849.1) | [MBG05719.1](https://www.ncbi.nlm.nih.gov/protein/MBG05719.1) | -       |
| Putative ID (16S)  | proteobacterium (91%) | Thalassospira GO-4 (86%) | - |
| Hit accession (16S)  | [AM259833.1](https://www.ncbi.nlm.nih.gov/nuccore/AM259833.1)    | [CP097807.1](https://www.ncbi.nlm.nih.gov/nuccore/CP097807.1) | -       |

### binning procedure ###
Unannotated scaffolds were binned.

```
#TODO
```

### annotation of genes with prodigal and barrnap ###
Using [prodigal](https://github.com/hyattpd/Prodigal) to identify CDS across all scaffolds for each bin:

```
~/prodigal-v2.6.3/prodigal -a Aphrocallistes_vastus_bacterial_Bin1.prot.fa -d Aphrocallistes_vastus_bacterial_Bin1.nucl.fa -f gff -o Aphrocallistes_vastus_bacterial_Bin1.gff -i Aphrocallistes_vastus_bacterial_Bin1.fasta
~/prodigal-v2.6.3/prodigal -a Aphrocallistes_vastus_bacterial_Bin2.prot.fa -d Aphrocallistes_vastus_bacterial_Bin2.nucl.fa -f gff -o Aphrocallistes_vastus_bacterial_Bin2.gff -i Aphrocallistes_vastus_bacterial_Bin2.fasta
```

The proteins were functionally annotated with the [blastKOALA](https://www.kegg.jp/blastkoala/) webserver. Pathway maps can be reconstructed with the [KEGG Mapper tool](https://www.kegg.jp/kegg/mapper/reconstruct.html) by uploading either of the `.user_ko_definition.txt` files.

Using [*E. coli* atpA](https://www.uniprot.org/uniprotkb/P0ABB0/entry) as the query, we identified the putative homologs of `atpA` and `atpB` in the two bacteria. These were then used to search against the `nr` database on [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to identify the closest relatives for assignment of the most precise names.

```
>ONT_TR_PE_scaff_693_297 hitby=sp|P0ABB0|ATPA_ECOLI e=0.0e+00 bits=687.6
MKQLSPAEISEVIKKKIGDINLDAEFRNEGTVLSVSDGIVRIHGLSAVQQGEMLSFSDGG
IGLALNLERDSVGAVVLGAYEHISEGGRVACTGKVLEVPVGEALIGRIVNPLGEAIDGKG
AIAATESSEPIEKIAPGVIWRQSVSQPLQTGLKAIDTMVPVGRGQRELIIGDRQTGKTAI
AVDTIINQKGQNCVCVYVAIGQKASSVAAVVRKLEEHGALEHTIIVSASASDSAALQYIA
PYSGCTMAEYFRDRGQDALIVYDDLTKQAWAYRQISLLLRRPPGREAYPGDVFYLHSRLL
ERASRVNPDYVEKQTNGEVKGKTGSLTALPIIETQAGDVTAFVPTNVISITDGQIFLETD
LFNAGIRPAMNAGLSVSRVGGAAQTKIIKKLGGGVRLALAQYRELAAFAQFASDLDDATR
KQLNIGKMVTELMKQKQYAPLKVSELALTLYAIERGFFAKVSVEKALKTEDDMQSAFAAK
HQALMDKIDAKPELTDEVDAELKAALTEFFENIG*
>Aphrocallistes_vastus_Bin2_222 hitby=sp|P0ABB0|ATPA_ECOLI e=0.0e+00 bits=559.7
MSLSGDDISALLKQNIADFDSIAEVSEIGHVISIGDGVARVYGLDLVRAGEMVEFESGVR
GMALNLETNTVGVVLFGDDRHVKEGETVKRTGSIVDLPVGKALLGRVVDALGNPIDGKGA
ITTEHRGLVDVKAPGIIPRQSVCEPVQTGIKALDTLVPIGRGQRELIIGDRQTGKSALAL
DTIINQRLANQSDDESKKLYCVYVVIGQKRSTVAQLARTLEEYGALEYTTIVAATASDPA
PLQFLAPYSGCAIGEYYRDNGMHALIVYDDLSKQANAYRQMSLLMRRPPGREAYPGDVFY
IHSRLLERAAKLSDERGGGSLTALPIIETQASDVSAYIPTNVISITDGQIFLQTDLFYKG
IRPAVNVGISVSRVGSAAQSKAMKSVAGSIKLELAQFREMEAFAQFASDLDSTTQNLLAR
GARLTELLKQKQYAPMPMEEQVVVIFAGVKGLLDTIPVERVTEFEESMLTRMRNSQSGVL
KDIRQSGKLSEENASKLTAALKRFAQEFLQGPAQDSTKGSSQGSIKGFGKAS*
```

Ribosomal RNA was identified using [barrnap](https://github.com/tseemann/barrnap). 16S was also identified in both bins, and used as a query on [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi). The 16S genes from both bins were also used as queries.

```
~/git/barrnap/bin/barrnap -o Aphrocallistes_vastus_bacterial_Bin1.barrnap.fa Aphrocallistes_vastus_bacterial_Bin1.fasta > Aphrocallistes_vastus_bacterial_Bin1.barrnap.gff
~/git/barrnap/bin/barrnap -o Aphrocallistes_vastus_bacterial_Bin2.barrnap.fa Aphrocallistes_vastus_bacterial_Bin2.fasta > Aphrocallistes_vastus_bacterial_Bin2.barrnap.gff
```



