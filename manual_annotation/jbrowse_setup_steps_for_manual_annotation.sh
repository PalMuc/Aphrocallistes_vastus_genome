#!/bin/bash

# APHROCALLISTES VASTUS PROCESSING STEPS

cd ~/genomes/aphrocallistes_vastus_PORI

~/samtools-1.9/samtools faidx Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted.fasta

# try to fix augustus format and jg/g bug
#~/git/genomeGTFtools/misc/augustus_to_gff3.py BRAKER2_ONT-RNA_minimap2_augustus.hints.gtf -x | sed s/jg/g/g > BRAKER2_ONT-RNA_minimap2_augustus.hints.gff3
# without jg to g
#~/git/genomeGTFtools/misc/augustus_to_gff3.py BRAKER2_ONT-RNA_minimap2_augustus.hints.gtf -x > BRAKER2_ONT-RNA_minimap2_augustus.hints.gff3
# jg to g is not the bug, filter two bad entries
~/git/genomeGTFtools/misc/augustus_to_gff3.py BRAKER2_ONT-RNA_minimap2_augustus.hints.gtf -x | grep -v "Parent=g15879" | grep -v "Parent=g1822" > BRAKER2_ONT-RNA_minimap2_augustus.hints.gff3

~/git/genomeGTFtools/misc/augustus_to_gff3.py BRAKER2_ONT-RNA_minimap2_braker.gtf -x | sed s/jg/g/g > BRAKER2_ONT-RNA_minimap2_braker.gff
~/git/genomeGTFtools/misc/augustus_to_gff3.py BRAKER2_PE-RNA_hisat2_augustus.hints.gtf -x | sed s/jg/g/g > BRAKER2_PE-RNA_hisat2_augustus.hints.gff


# fix stringtie format
~/git/genomeGTFtools/misc/stringtie_gtf_to_gff3.py STRG_PE-RNA_hisat2_stranded.gtf > STRG_PE-RNA_hisat2_stranded.gff
~/git/genomeGTFtools/misc/stringtie_gtf_to_gff3.py STRG_ONT-RNA_minimap2.gtf > STRG_ONT-RNA_minimap2.gff



# add raw ONT reads
~/samtools-1.9/samtools index Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.sorted.bam


# blast vs models

~/diamond-linux64/diamond blastp -q BRAKER2_ONT-RNA_augustus.hints_proteins.fasta -d ~/db/model_organism_uniprot.fasta -o BRAKER2_ONT-RNA_augustus.hints_proteins.vs_models.tab

~/git/genomeGTFtools/blast2genomegff.py -b BRAKER2_ONT-RNA_augustus.hints_proteins.vs_models.tab -d ~/db/model_organism_uniprot.fasta -p blastp -g BRAKER2_ONT-RNA_minimap2_augustus.hints.gff3 --add-accession --add-description -P -S -x -K -G -F "." > BRAKER2_ONT-RNA_augustus.hints_proteins.vs_models.gff

# blast vs demosponges
~/diamond-linux64/diamond blastp -q BRAKER2_ONT-RNA_augustus.hints_proteins.fasta -d ~/est/porifera/demo_genomes_combined.fasta -o BRAKER2_ONT-RNA_augustus.hints_proteins.vs_demos.tab
~/git/genomeGTFtools/blast2genomegff.py -b BRAKER2_ONT-RNA_augustus.hints_proteins.vs_demos.tab -d ~/est/porifera/demo_genomes_combined.fasta -p blastp -g BRAKER2_ONT-RNA_minimap2_augustus.hints.gff3 -P -x -K -G -F "." -M 5 > BRAKER2_ONT-RNA_augustus.hints_proteins.vs_demos.gff

# blast vs glass sponges

cat hyalonema_populiferum_prots.fasta sympagella_nux_prots.fasta rosella_fibulata_prots.fasta > glass_sponge_combined.fasta
~/diamond-linux64/diamond makedb --in glass_sponge_combined.fasta --db glass_sponge_combined.fasta
~/diamond-linux64/diamond blastp -q BRAKER2_ONT-RNA_augustus.hints_proteins.fasta -d ~/est/porifera/glass_sponge_combined.fasta -o BRAKER2_ONT-RNA_augustus.hints_proteins.vs_glass.tab 
~/git/genomeGTFtools/blast2genomegff.py -b BRAKER2_ONT-RNA_augustus.hints_proteins.vs_glass.tab -d ~/est/porifera/glass_sponge_combined.fasta -p blastp -g BRAKER2_ONT-RNA_minimap2_augustus.hints.gff3 -P -x -K -G -F "." > BRAKER2_ONT-RNA_augustus.hints_proteins.vs_glass.gff
~/git/genomeGTFtools/blast2genomegff.py -b BRAKER2_ONT-RNA_augustus.hints_proteins.vs_glass.tab -d ~/est/porifera/glass_sponge_combined.fasta -p blastp -g BRAKER2_ONT-RNA_minimap2_augustus.hints.gff3 -P -x -K -G -F "." -M 5 > BRAKER2_ONT-RNA_augustus.hints_proteins.vs_glass.gff


# check for trans splice
#>STRG_ONT-RNA_minimap2_transcripts.533.1.exon2
#TCATTTGATTTACTCAGCATTTATTCTATTCAATTTCAAGCATTAGAAAGTTTATTCTAACTTCAGATTACTGTTCTTACCTTAGATTTTCAG
#>STRG_ONT-RNA_minimap2_transcripts.535.1.exon2
#AGCTTCAATATATTCTTATGCTGGCTCAGTTCCTATATAAATAGCATAAATTGATGTGTCTAAACTCAATCGATTTAATCAGCATTTGTTCTATTCAATTTCAAGCTTTCAAAAGTTTATTCTAACTTCATATTACTGTTCTTATCTCAGATTTTCAG
#>STRG_ONT-RNA_minimap2_transcripts.546.1.exon2
#AATCGATTTATTTGTTCTAAATTTGTTCTAAACTCTTTATTTGTTCTAAACTCAATCGATTCATTTGTTCTATTCATTTATCGAGCTTTCAAAAGTTTAATCTAAGCTTCAGATTACTGTTCTTATCTTTGATTTTCAG


# manual editing tracks
cat BRAKER2_PE-RNA_augustus.hints.unique.gff3 | grep -v "Parent=g13021" | grep -v "Parent=g6171" | grep -v "Parent=g8452" | grep -v "Parent=g5955" > BRAKER2_PE-RNA_augustus.hints.unique.gff3.temp
mv BRAKER2_PE-RNA_augustus.hints.unique.gff3.temp BRAKER2_PE-RNA_augustus.hints.unique.gff3



~/git/genomeGTFtools/misc/stringtie_gtf_to_gff3.py STRG_ONT-RNA_minimap2_f-0.3.unique.gtf | sed s/STRG_//g | sed s/_transcripts/_tx/g > STRG_ONT-RNA_minimap2_f-0.3.unique.gff

~/git/genomeGTFtools/misc/stringtie_gtf_to_gff3.py STRG_PE-RNA_hisat2_non-stranded_f-0.3.unique.gtf | sed s/STRG_//g | sed s/non-stranded_non-stranded_transcripts/NS_tx/g > STRG_PE-RNA_hisat2_non-stranded_f-0.3.unique.gff

~/git/genomeGTFtools/misc/stringtie_gtf_to_gff3.py STRG_PE-RNA_hisat2_stranded_f-0.3.unique.gtf | sed s/STRG_//g | sed s/stranded_transcripts/Sr_tx/g > STRG_PE-RNA_hisat2_stranded_f-0.3.unique.gff


~/minimap2-2.14_x64-linux/minimap2 -a -x splice --secondary=no Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted.fasta Avastus_individual-1_trinity.fasta | ~/samtools-1.8/samtools sort - -o Avastus_individual-1_trinity.bam
~/samtools-1.8/samtools index Avastus_individual-1_trinity.bam



# blast pinfish tx

~/gffread-0.12.1.Linux_x86_64/gffread PINFISH_pipeline_raw-reads_clustered_transcripts.unique.gff -g ../Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted.fasta -w PINFISH_pipeline_raw-reads_clustered_transcripts.unique.fasta
~/diamond-linux64/diamond blastx -q PINFISH_pipeline_raw-reads_clustered_transcripts.unique.fasta -d ~/est/porifera/glass_sponge_combined.fasta -o PINFISH_pipeline_raw-reads_clustered_transcripts.unique.vs_glass.tab
~/git/genomeGTFtools/blast2genomegff.py -b PINFISH_pipeline_raw-reads_clustered_transcripts.unique.vs_glass.tab -d ~/est/porifera/glass_sponge_combined.fasta -p blastx -g PINFISH_pipeline_raw-reads_clustered_transcripts.unique.gff -P -M 3 > PINFISH_pipeline_raw-reads_clustered_transcripts.unique.vs_glass.gff


#
# CREATE JBROWSE
#

cd ~/genomes/jbrowse_data/aphrocallistes
ln -s /mnt/data/genomes/aphrocallistes_vastus/Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted.fasta
ln -s /mnt/data/genomes/aphrocallistes_vastus/Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted.fasta.fai

/var/www/html/jbrowse/bin/prepare-refseqs.pl --indexed_fasta Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted.fasta --out ./

ln -s /mnt/data/genomes/aphrocallistes_vastus/BRAKER2_ONT-RNA_minimap2_augustus.hints.gff3

/var/www/html/jbrowse/bin/flatfile-to-json.pl --gff BRAKER2_ONT-RNA_minimap2_augustus.hints.gff3 --trackType CanvasFeatures  --trackLabel BRAKER2_ONT-RNA_minimap2_augustus --out ./

ln -s /mnt/data/genomes/aphrocallistes_vastus/STRG_PE-RNA_hisat2_stranded.gff
ln -s /mnt/data/genomes/aphrocallistes_vastus/STRG_ONT-RNA_minimap2.gff
/var/www/html/jbrowse/bin/flatfile-to-json.pl --gff STRG_PE-RNA_hisat2_stranded.gff --trackType CanvasFeatures  --trackLabel STRG_PE-RNA_hisat2_stranded --out ./
/var/www/html/jbrowse/bin/flatfile-to-json.pl --gff STRG_ONT-RNA_minimap2.gff --trackType CanvasFeatures  --trackLabel STRG_ONT-RNA_minimap2 --out ./

ln -s /mnt/data/genomes/aphrocallistes_vastus/Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.sorted.bam
ln -s /mnt/data/genomes/aphrocallistes_vastus/Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.sorted.bam.bai 

ln -s /mnt/data/genomes/aphrocallistes_vastus/Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.PE-RNA_hisat2_stranded.fwd.bg.cov.bw 
ln -s /mnt/data/genomes/aphrocallistes_vastus/Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.PE-RNA_hisat2_stranded.rev.bg.cov.bw 


ln -s /mnt/data/genomes/aphrocallistes_vastus/BRAKER2_ONT-RNA_BRAKER2_ONT-RNA_protein_alignment_gth.gff3
/var/www/html/jbrowse/bin/flatfile-to-json.pl --out ./ --gff BRAKER2_ONT-RNA_BRAKER2_ONT-RNA_protein_alignment_gth.gff3 --trackType CanvasFeatures --trackLabel ONT-RNA_protein_alignment_gth

# add blast
ln -s /mnt/data/genomes/aphrocallistes_vastus/BRAKER2_ONT-RNA_augustus.hints_proteins.vs_models.gff

/var/www/html/jbrowse/bin/flatfile-to-json.pl --out ./ --gff BRAKER2_ONT-RNA_augustus.hints_proteins.vs_models.gff --trackType CanvasFeatures --trackLabel BRAKER2_ONT_augustus_vs_models

ln -s /mnt/data/genomes/aphrocallistes_vastus/BRAKER2_ONT-RNA_augustus.hints_proteins.vs_glass.gff
/var/www/html/jbrowse/bin/flatfile-to-json.pl --out ./ --gff BRAKER2_ONT-RNA_augustus.hints_proteins.vs_glass.gff --trackType CanvasFeatures --trackLabel BRAKER2_ONT_augustus_vs_glass_tx

ln -s /mnt/data/genomes/aphrocallistes_vastus/BRAKER2_ONT-RNA_augustus.hints_proteins.vs_demos.gff
/var/www/html/jbrowse/bin/flatfile-to-json.pl --out ./ --gff BRAKER2_ONT-RNA_augustus.hints_proteins.vs_demos.gff --trackType CanvasFeatures --trackLabel BRAKER2_ONT_augustus_vs_demo


###
### RE CREATE JBROWSE with UNIQUE NAMES
###

ln -s /mnt/data/genomes/aphrocallistes_vastus/manual_curation_tracks/BRAKER2_ONT-RNA_minimap2_augustus.unique.gff3 

ln -s /mnt/data/genomes/aphrocallistes_vastus/manual_curation_tracks/BRAKER2_PE-RNA_augustus.hints.unique.gff3

/var/www/html/jbrowse/bin/flatfile-to-json.pl --gff BRAKER2_ONT-RNA_minimap2_augustus.unique.gff3 --trackType CanvasFeatures  --trackLabel BRAKER2_ONT-RNA_augustus --out ./

/var/www/html/jbrowse/bin/flatfile-to-json.pl --gff BRAKER2_PE-RNA_augustus.hints.unique.gff3 --trackType CanvasFeatures  --trackLabel BRAKER2_PE-RNA_augustus --out ./


ln -s /mnt/data/genomes/aphrocallistes_vastus/manual_curation_tracks/STRG_ONT-RNA_minimap2_f-0.3.unique.gff
/var/www/html/jbrowse/bin/flatfile-to-json.pl --gff STRG_ONT-RNA_minimap2_f-0.3.unique.gff --trackType CanvasFeatures  --trackLabel STRG_ONT-RNA_f-0.3 --out ./

ln -s /mnt/data/genomes/aphrocallistes_vastus/manual_curation_tracks/STRG_PE-RNA_hisat2_non-stranded_f-0.3.unique.gff
/var/www/html/jbrowse/bin/flatfile-to-json.pl --gff STRG_PE-RNA_hisat2_non-stranded_f-0.3.unique.gff --trackType CanvasFeatures  --trackLabel STRG_PE-RNA_NS_f-0.3 --out ./

ln -s /mnt/data/genomes/aphrocallistes_vastus/manual_curation_tracks/STRG_PE-RNA_hisat2_stranded_f-0.3.unique.gff
/var/www/html/jbrowse/bin/flatfile-to-json.pl --gff STRG_PE-RNA_hisat2_stranded_f-0.3.unique.gff --trackType CanvasFeatures  --trackLabel STRG_PE-RNA_Sr_f-0.3 --out ./

ln -s /mnt/data/genomes/aphrocallistes_vastus/manual_curation_tracks/PINFISH_pipeline_raw-reads_clustered_transcripts.unique.gff
/var/www/html/jbrowse/bin/flatfile-to-json.pl --gff PINFISH_pipeline_raw-reads_clustered_transcripts.unique.gff --trackType CanvasFeatures  --trackLabel pinfish-raw-clustered --out ./









