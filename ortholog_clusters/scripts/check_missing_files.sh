#/bin/bash

ls clusters_holozoa_allprots_v10/homologs/*fasta > clusters_holozoa_allprots_v10_homologs.txt
sed -i 's/clusters_holozoa_allprots_v10\/homologs\///' clusters_holozoa_allprots_v10_homologs.txt
ls clusters_holozoa_allprots_v10/paralogs/*fasta > clusters_holozoa_allprots_v10_paralogs.txt
sed -i 's/clusters_holozoa_allprots_v10\/paralogs\///' clusters_holozoa_allprots_v10_paralogs.txt
cat clusters*txt > clusters_holozoa_allprots_v10_all.txt


ls mafft/paralogs/*aln > aln_list.paralogs.txt
sed -i 's/mafft\/paralogs\///' aln_list.paralogs.txt
sed -i 's/\.aln//' aln_list.paralogs.txt 

ls mafft/homologs/*aln > aln_list.homologs.txt
sed -i 's/mafft\/homologs\///' aln_list.homologs.txt
sed -i 's/\.aln//' aln_list.homologs.txt 

ls fasttree/paralogs/*tre > tre_list.paralogs.txt
sed -i 's/fasttree\/paralogs\///' tre_list.paralogs.txt 
sed -i 's/\.aln\.tre//' tre_list.paralogs.txt 

ls fasttree/homologs/*tre > tre_list.homologs.txt
sed -i 's/fasttree\/homologs\///' tre_list.homologs.txt 
sed -i 's/\.aln\.tre//' tre_list.homologs.txt 

ls pfam/paralogs/*tab > pfam_list.paralogs.txt
sed -i 's/pfam\/paralogs\///' pfam_list.paralogs.txt
sed -i 's/\.pfam.tab/.fasta/' pfam_list.paralogs.txt 

ls pfam/homologs/*tab > pfam_list.homologs.txt
sed -i 's/pfam\/homologs\///' pfam_list.homologs.txt
sed -i 's/\.pfam.tab/.fasta/' pfam_list.homologs.txt 

cat tre_list*txt > tre_list.all.txt
cat aln_list*txt > aln_list.all.txt
cat pfam_list*txt > pfam_list.all.txt




grep -Ff aln_list.all.txt -v clusters_holozoa_allprots_v10_all.txt > missing_alignments.txt
grep -Ff tre_list.all.txt -v clusters_holozoa_allprots_v10_all.txt > missing_trees.txt
grep -Ff pfam_list.all.txt -v clusters_holozoa_allprots_v10_all.txt > missing_pfams.txt


