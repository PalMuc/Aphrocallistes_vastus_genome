# steps to make final gene annotation from manual annotation #

### starting from BRAKER annotation ###

`~/git/genomeGTFtools/misc/augustus_to_gff3.py BRAKER2_ONT-RNA_minimap2_augustus.hints.gtf -x | sed s/jg/g/g > BRAKER2_ONT-RNA_minimap2_augustus.hints.gff3`

`./prepare_names_for_spreadsheet.py BRAKER2_ONT-RNA_minimap2_augustus.hints.gff3 > names_for_spreadsheet`

### convert annotation table into rough GFF ###

```
cd ~/genomes/aphrocallistes_vastus_PORI/

mv ~/Downloads/Avas\ manual\ annotation\ -\ annotation.tsv Avas_manual_annotation_2022-02-11.tsv

./clean_avas_annotation_table.py Avas_manual_annotation_2022-02-11.tsv > Avas_manual_annotation_2022-02-11.fix.tsv
# removing comment lines, header line, was_checked column, and any columns after comment_for_annotators
# counted 20644 lines
# annotations by:
# m : 11433
# w : 3583
# r : 3280
# n : 717
# s : 552
# k : 494
# c : 463
# e : 114
#  : 9

# from WRF https://bitbucket.org/wrf/genome-reannotations/src/master/
```

`~/git/genome-reannotations/remake_genome_annotation.py -m manually_edited_genes_revised.gff -i Avas_manual_annotation_2022-02-11.fix.tsv -g BRAKER2_ONT-RNA.revised.sort.gff3 -a BRAKER2_ref_proteins_aln.gff3 BRAKER2_PE-RNA.gff3 PINFISH_stranded.gff PINFISH_stepwise.gff PINFISH_pipeline.gff STRG_PE-RNA_stranded.gff STRG_PE-RNA_hisat2.gff STRG_ONT-RNA.gff STRG_ONT-RNA_minimap2.gff -s Avas.v1.29 -f Avas.v1.29_names.tab --clean-scaffold-text Aphrocalllistes_vastus_HiC-scaffold_ --replace-scaffold-text Avas.s > Avas.v1.29_annotations.gff 2> Avas.v1.29_annotations.log`

`fix_avas_rrna_gff.py Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_rRNA_revised.gff`

### transcript translate and add manual nucleotides ###

```
# from https://github.com/gpertea/gffread
~/gffread-0.12.7.Linux_x86_64/gffread Avas.v1.29_annotations.gff -g Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted.fasta -w Avas.v1.29_annotations.nucl.fasta -E

#$  grep ">" Avas.v1.29_annotations.nucl.fasta -c 
#20052
#$  grep "CDS=" Avas.v1.29_annotations.nucl.fasta -c 
#17928
#$  grep "CDS=1" Avas.v1.29_annotations.nucl.fasta -c
#17764

./replace_manual_cds_nucls.py Avas.v1.29_annotations.gff manual_genes_to_map.fasta Avas.v1.29_annotations.nucl.fasta > Avas.v1.29_annotations.corr_nucl.fasta

prottrans.py -r -v -x Avas.v1.29_annotations.corr_nucl.fasta > Avas.v1.29_annotations.corr_prot.fasta
prottrans.py -r -v -x -t f Avas.v1.29_annotations.corr_nucl.fasta > Avas.v1.29_annotations.full_frame_prot.fasta

./fix_final_avas_annotation.py Avas.v1.29_annotations.full_frame_prot.fasta Avas.v1.29_annotations.corr_prot.fasta Avas.v1.29_annotations.gff > Avas.v1.29_annotations.corr.gff

# diagnostic
~/git/genome-reannotations/remake_genome_annotation.py -g Avas.v1.29_annotations.corr.gff -i null
# Getting main gene annotation from -g
# Reading GFF file Avas.v1.29_annotations.corr.gff
# File contained:
#exon	50946
#CDS	49636
#mRNA	20038
#gene	19578
#rRNA	22
# Collected exons for 20060 transcripts
```

### tidy up datasets ###

```
cat Avas.v1.29_annotations.corr.gff | sed s/Aphrocalllistes/Aphrocallistes/g > v1.29/Avas.v1.29_annotations.fr.gff
cat v1.29/Avas.v1.29_annotations.fr.gff | gzip > v1.29/Avas.v1.29_annotations.fr.gff.gz
cat Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted.fasta | sed s/Aphrocalllistes/Aphrocallistes/g > v1.29/Avas.v1.29_scaffolds.fasta
cat v1.29/Avas.v1.29_scaffolds.fasta | gzip > v1.29/Avas.v1.29_scaffolds.fasta.gz
cp Avas.v1.29_annotations.corr_nucl.fasta v1.29/Avas.v1.29_annotations.transcripts.fr.fasta
cat v1.29/Avas.v1.29_annotations.transcripts.fr.fasta | gzip > v1.29/Avas.v1.29_annotations.transcripts.fr.fasta.gz
prottrans.py -r -v -x Avas.v1.29_annotations.transcripts.fr.fasta > v1.29/Avas.v1.29_annotations.prot.fr.fasta
cat v1.29/Avas.v1.29_annotations.prot.fr.fasta | gzip > v1.29/Avas.v1.29_annotations.prot.fr.fasta.gz

# remove old version, and relink to new version
rm Avas.current_v.annotations.gff ; ln -s v1.29/Avas.v1.29_annotations.fr.gff Avas.current_v.annotations.gff
```

### remake some tracks for browser ###

```
# some functional annotation
~/diamond-v2.0.13/diamond-v2.0.13 blastx -q v1.29/Avas.v1.29_annotations.transcripts.fr.fasta -d ~/db/model_organism_uniprot.w_cnido.fasta -o v1.29/Avas.v1.29_annotations.blast_vs_models.tab
~/git/genomeGTFtools/blast2genomegff.py -b v1.29/Avas.v1.29_annotations.blast_vs_models.tab -d ~/db/model_organism_uniprot.w_cnido.fasta -g v1.29/Avas.v1.29_annotations.fr.gff -S --add-description --add-accession > v1.29/Avas.v1.29_annotations.blast_vs_models.gff

~/samtools-1.14/samtools fastq -n -s Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.fastq Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.sorted.bam > Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.fastq
~/minimap2-2.23_x64-linux/minimap2 -a -x splice --secondary=no v1.29/Avas.v1.29_scaffolds.fasta Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.fastq | ~/samtools-1.14/samtools sort - -o v1.29/Avas.v1.29_scaffolds.ONT-RNA_minimap2.bam
~/samtools-1.14/samtools index Avas.v1.29_scaffolds.ONT-RNA_minimap2.bam
```

### comparison among annotations ###
Using BUSCO 5.4.2 for dataset `metazoa_odb10` with 954 proteins:

| annotation | complete | single | duplicated | fragmented | missing |
| :--- | --- | --- | --- | --- | --- |
| Avas.v1.29_annotations | 80.4% | 78.5% | 1.9% | 2.3% | 17.3% |
| BRAKER2_ONT-RNA_augustus.hints_proteins | 77.8% | 74.7% | 3.1% | 3.4% | 18.8% |
| PINFISH_pipeline | 69.9% | 54.6% | 15.3% | 4.7% | 25.4% |
| STRG_ONT-RNA_minimap2_transdecoder| 76.6% | 58.5% | 18.1% | 3.1% | 20.3% |


