#!/bin/bash


cat 02* 03* 04* 05* 06* 07* 13* 14* > A_choano_pori_allprots_v7.fasta.gz
cat 02* 04* 05* 06* 07* 13* 14* > B_mbre_sponges_allprots_v7.fasta.gz
cat 02* 03* 04* 05* 06* 07* 08* 09* 10* 11* 12* 13* 14* > C_holozoa_allprots_v7.fasta.gz
cat 02* 03* 04* 05* 06* 07* 08* 09* 10* 11* 12* HEX06* HEX07* HEX10* HEX12*  > D_holozoa_allprots_v7_reduced.fasta.gz
cat 02* 04* 05* 06* 07* 08* 09* 10* 11* 12* HEX06* HEX07*  > E_holozoa_allprots_v8.fasta.gz

