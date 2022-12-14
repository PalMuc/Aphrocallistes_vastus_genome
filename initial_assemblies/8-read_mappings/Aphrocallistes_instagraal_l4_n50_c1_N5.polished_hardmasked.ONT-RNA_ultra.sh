TESTDIR=/dss/dssfs02/lwp-dss-0001/uk213/uk213-dss-0000/di29vos/tools/anaconda3/envs/genomics/ultra/test


uLTRA pipeline $TESTDIR/SIRV_genes_C_170612a.gtf  \
               $TESTDIR/SIRV_genes.fasta  \
               $TESTDIR/reads.fa outfolder/ --t 8
               
               
uLTRA align  genome.fasta  reads.[fa/fq] outfolder/  --ont --t 48   # ONT cDNA reads using 48 cores
uLTRA align  genome.fasta  reads.[fa/fq] outfolder/  --isoseq --t 48 # PacBio isoseq reads
uLTRA align  genome.fasta  reads.[fa/fq] outfolder/  --k 14  --t 48 # PacBio dRNA reads or reads with >10-12% error rate
