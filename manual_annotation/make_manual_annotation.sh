# steps to make final gene annotation from manual annotation #

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
~/git/genome-reannotations/remake_genome_annotation.py -m manually_edited_genes_revised.gff -i Avas_manual_annotation_2022-02-11.fix.tsv -g BRAKER2_ONT-RNA.revised.sort.gff3 -a BRAKER2_ref_proteins_aln.gff3 BRAKER2_PE-RNA.gff3 PINFISH_stranded.gff PINFISH_stepwise.gff PINFISH_pipeline.gff STRG_PE-RNA_stranded.gff STRG_PE-RNA_hisat2.gff STRG_ONT-RNA.gff STRG_ONT-RNA_minimap2.gff -s Avas.v1.29 -f Avas.v1.29_names.tab --clean-scaffold-text Aphrocalllistes_vastus_HiC-scaffold_ --replace-scaffold-text Avas.s > Avas.v1.29_annotations.gff 2> Avas.v1.29_annotations.log

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

