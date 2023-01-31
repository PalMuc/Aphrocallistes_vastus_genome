# phylogeny of hexactinellids #
Most of the data is from the 2017 RV SONNE cruise as part of Project PoribacNewZ. The 454 transcriptome of [Euplectella aspergillum](https://marinegenomics.oist.jp/kairou/viewer/info?project_id=62) is from the [OIST Marine Genomics Unit](https://marinegenomics.oist.jp/gallery/gallery/index), and the three other transcriptomes come from [Whelan et al 2015](https://figshare.com/articles/dataset/Error_signal_and_the_placement_of_Ctenophora_sister_to_all_other_animals/1334306).

## generate supermatrix ##

Most of the code can be found in [WRF's supermatrix repo](https://github.com/wrf/supermatrix).

Make overall supermatrix:

```
./collate_busco_result_fasta.py -s species_to_dir_list.tab -o combined_buscos/

for FILE in combined_buscos/*.faa ; do mafft $FILE > $FILE.aln ; done

~/git/supermatrix/join_alignments.py -a combined_buscos/*.aln -d "|" -u hexa_sponge_combined_buscos.aln
# Finished parsing 926 alignments for 22 total taxa
```

Sort and make plot for all taxa for all genes:

```
~/git/supermatrix/reorder_matrix_by_cov.py -a hexa_sponge_combined_buscos.aln -p hexa_sponge_combined_buscos.aln.partition.txt -o hexa_sponge_combined_buscos.sorted.aln
# alignment has 926 partitions with 526348 sites

~/git/supermatrix/check_supermatrix_alignments.py -a hexa_sponge_combined_buscos.sorted.aln -p hexa_sponge_combined_buscos.sorted.aln.partition.txt -m hexa_sponge_combined_buscos.occ_matrix.tab --percent
# Matrix has 101 (0.50%) complete and 9928 (48.73%) empty out of 20372 total genes
Rscript ~/git/supermatrix/draw_matrix_occupancy.R hexa_sponge_combined_buscos.occ_matrix.tab
```

Filter for 50% coverage and require certain low coverage taxa:

```
~/git/supermatrix/filter_supermatrix.py -a hexa_sponge_combined_buscos.aln -p hexa_sponge_combined_buscos.aln.partition.txt -r Hyalonema_populiferum Walteria_leuckarti Corbitella_sp Farrea_similaris -o hexa_sponge_combined_buscos.filt.aln
~/git/supermatrix/reorder_matrix_by_cov.py -a hexa_sponge_combined_buscos.filt.aln -p hexa_sponge_combined_buscos.filt.aln.partition.txt -o hexa_sponge_combined_buscos.filt.sorted.aln
~/git/supermatrix/check_supermatrix_alignments.py -a hexa_sponge_combined_buscos.filt.sorted.aln -p hexa_sponge_combined_buscos.filt.sorted.aln.partition.txt -m hexa_sponge_combined_buscos.filt.occ_matrix.tab --percent
Rscript ~/git/supermatrix/draw_matrix_occupancy.R hexa_sponge_combined_buscos.filt.occ_matrix.tab
```

Make draft tree:

`FastTree hexa_sponge_combined_buscos.filt.sorted.aln > hexa_sponge_combined_buscos.filt.sorted.tree`

Run IQTree:

`iqtree2 -s hexa_sponge_combined_buscos.filt.sorted.aln -mset LG -b 100 -nt 32`

## process the BUSCO dataset ##
for various analyses. The [human NCBI ID](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9606) is `9606`, [mouse NCBI ID](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10090) is `10090`:

```
zcat refseq_db.faa.gz | grep ">" > refseq_db.all_names
grep _9606_ refseq_db.all_names > refseq_db.human_names
getAinB.py refseq_db.human_names refseq_db.faa.gz > refseq_db.human.faa
blastp -query refseq_db.human.faa -db ../human_uniprot.fasta -outfmt 6 -max_target_seqs 1 > refseq_db.human.vs_uniprot.tab

grep _10090_ refseq_db.all_names > refseq_db.mouse_names
getAinB.py refseq_db.mouse_names refseq_db.faa.gz > refseq_db.mouse.faa
getAinB.py _MOUSE uniprot_sprot.fasta -s > mouse_uniprot.fasta
blastp -query refseq_db.mouse.faa -db ../mouse_uniprot.fasta -outfmt 6 -max_target_seqs 1 > refseq_db.mouse.vs_uniprot.tab
```

Pull list of BUSCOs not found in any species:

```
./combine_busco_table_results.py -s species_to_dir_list.tab > hexa_species_all_busco_results.tab
# Wrote 954 genes for 22 species

grep -E -f pori_missing_buscos.txt ~/db/metazoa_odb10/refseq_db.human.vs_uniprot.tab
522988at33208_9606_0:0010df	sp|Q96IV0|NGLY1_HUMAN	100.000	654	0	0	1	654	1	654	0.0	1357
344346at33208_9606_0:0015f1	sp|Q92552|RT27_HUMAN	96.729	428	0	1	1	428	1	414	0.0	843
608669at33208_9606_0:00266e	sp|O95900|TRUB2_HUMAN	100.000	331	0	0	1	331	1	331	0.0	681
346541at33208_9606_0:003067	sp|O00255|MEN1_HUMAN	93.455	657	1	1	1	657	1	615	0.0	1241
257141at33208_9606_0:0038f4	sp|O94923|GLCE_HUMAN	100.000	617	0	0	7	623	1	617	0.0	1292
558385at33208_9606_0:003a9e	sp|Q9Y2V0|CO041_HUMAN	97.909	287	0	1	1	287	1	281	0.0	580
556039at33208_9606_0:003c27	sp|O96000|NDUBA_HUMAN	100.000	172	0	0	1	172	1	172	1.56e-127	354
359532at33208_9606_0:003c62	sp|Q8TBB5|KLDC4_HUMAN	98.857	525	1	1	1	525	1	520	0.0	1067
643954at33208_9606_0:004cf6	sp|Q9BYN8|RT26_HUMAN	100.000	205	0	0	1	205	1	205	4.26e-147	406

grep -E -f pori_missing_buscos.txt ~/db/metazoa_odb10/refseq_db.mouse.vs_uniprot.tab
558385at33208_10090_0:0007cc	sp|Q3U4G0|CO041_MOUSE	99.174	242	2	0	1	242	1	242	3.18e-180	500
643954at33208_10090_0:000976	sp|Q80ZS3|RT26_MOUSE	100.000	200	0	0	1	200	1	200	2.05e-142	394
608669at33208_10090_0:000b89	sp|Q91WG3|TRUB2_MOUSE	100.000	331	0	0	1	331	1	331	0.0	680
361216at33208_10090_0:001122	sp|Q6PDJ6|FBX42_MOUSE	100.000	717	0	0	1	717	1	717	0.0	1470
526839at33208_10090_0:0013c7	sp|Q921S7|RM37_MOUSE	100.000	423	0	0	1	423	1	423	0.0	872
495802at33208_10090_0:001d18	sp|Q9CPY1|RM51_MOUSE	100.000	128	0	0	1	128	1	128	1.72e-92	262
359532at33208_10090_0:0027ae	sp|Q921I2|KLDC4_MOUSE	99.315	584	4	0	1	584	1	584	0.0	1196
344346at33208_10090_0:003e26	sp|Q8BK72|RT27_MOUSE	91.009	456	0	1	1	456	1	415	0.0	831
557694at33208_10090_0:003f65	sp|Q8R3J4|MTEF3_MOUSE	100.000	412	0	0	1	412	1	412	0.0	834
239149at33208_10090_0:003fca	sp|Q9D0G0|RT30_MOUSE	100.000	442	0	0	1	442	1	442	0.0	913
522988at33208_10090_0:0040f7	sp|Q9JI78|NGLY1_MOUSE	100.000	651	0	0	1	651	1	651	0.0	1343
560566at33208_10090_0:004dd2	sp|P19426|NELFE_MOUSE	97.867	375	0	1	1	367	1	375	0.0	730
556039at33208_10090_0:004df7	sp|Q9DCS9|NDUBA_MOUSE	100.000	176	0	0	1	176	1	176	9.88e-132	365
640187at33208_10090_0:005253	sp|Q9JHK4|PGTA_MOUSE	28.070	171	88	5	89	250	50	194	1.10e-06	50.4

grep -E -f hexa_missing_buscos.txt ~/db/metazoa_odb10/refseq_db.human.vs_uniprot.tab
40200at33208_9606_0:000a8a	sp|Q9UBZ9|REV1_HUMAN	99.920	1251	0	1	109	1358	1	1251	0.0	2603
601875at33208_9606_0:001395	sp|Q4G0J3|LARP7_HUMAN	97.815	595	0	1	1	595	1	582	0.0	1172
502584at33208_9606_0:001479	sp|Q16585|SGCB_HUMAN	100.000	318	0	0	1	318	1	318	0.0	659
334272at33208_9606_0:0017aa	sp|Q8N302|AGGF1_HUMAN	100.000	714	0	0	1	714	1	714	0.0	1484
590347at33208_9606_0:0018d8	sp|Q0VDI3|TM267_HUMAN	100.000	215	0	0	1	215	1	215	4.13e-159	437
360379at33208_9606_0:00206e	sp|Q1RMZ1|BMT2_HUMAN	100.000	405	0	0	1	405	1	405	0.0	845
393219at33208_9606_0:002491	sp|Q9BTY7|HGH1_HUMAN	100.000	390	0	0	1	390	1	390	0.0	760
625091at33208_9606_0:0024b3	sp|Q9Y2Q9|RT28_HUMAN	100.000	187	0	0	1	187	1	187	3.61e-138	382
549988at33208_9606_0:0024fb	sp|Q96E11|RRFM_HUMAN	100.000	262	0	0	22	283	1	262	0.0	537
412930at33208_9606_0:002a01	sp|Q96KC8|DNJC1_HUMAN	97.105	380	8	1	1	377	1	380	0.0	757
412930at33208_9606_0:002a01	sp|Q96KC8|DNJC1_HUMAN	29.630	54	37	1	305	358	473	525	8.2	28.5
493505at33208_9606_0:002f69	sp|Q96GX9|MTNB_HUMAN	98.253	229	4	0	31	259	14	242	1.17e-171	473
203698at33208_9606_0:003223	sp|Q8IY37|DHX37_HUMAN	100.000	1129	0	0	1	1129	1	1129	0.0	2299
641073at33208_9606_0:00350a	sp|Q9UHK0|NUFP1_HUMAN	100.000	495	0	0	1	495	1	495	0.0	1024
449056at33208_9606_0:0038ad	sp|P51948|MAT1_HUMAN	97.909	287	6	0	8	294	23	309	0.0	571
613522at33208_9606_0:003b3f	sp|Q9H2W6|RM46_HUMAN	100.000	279	0	0	1	279	1	279	0.0	570
524148at33208_9606_0:003bae	sp|Q9BQ70|TCF25_HUMAN	99.174	605	0	1	1	605	1	600	0.0	1224
648300at33208_9606_0:0042f3	sp|Q96MW1|CCD43_HUMAN	100.000	143	0	0	1	143	1	143	8.28e-100	285
495329at33208_9606_0:004558	sp|Q9NQT4|EXOS5_HUMAN	100.000	235	0	0	1	235	1	235	5.33e-175	479
534172at33208_9606_0:00455f	sp|Q9BYD3|RM04_HUMAN	100.000	246	0	0	1	246	1	246	9.96e-179	497
619250at33208_9606_0:00495b	sp|Q5BKX5|CS054_HUMAN	100.000	325	0	0	1	325	1	325	0.0	647
399439at33208_9606_0:004d16	sp|P57081|WDR4_HUMAN	100.000	412	0	0	1	412	1	412	0.0	845

grep -E -f hexa_missing_buscos.txt ~/db/metazoa_odb10/refseq_db.mouse.vs_uniprot.tab
40200at33208_10090_0:00008f	sp|Q920Q2|REV1_MOUSE	100.000	1249	0	0	58	1306	1	1249	0.0	2598
503828at33208_10090_0:00029d	sp|Q8R2X8|GO45_MOUSE	98.293	410	0	1	1	410	1	403	0.0	826
279992at33208_10090_0:0003b8	sp|Q8C735|LIN9_MOUSE	100.000	542	0	0	17	558	1	542	0.0	1130
360753at33208_10090_0:00053c	sp|A2APY7|NDUF5_MOUSE	100.000	343	0	0	1	343	1	343	0.0	715
549988at33208_10090_0:00081b	sp|Q9D6S7|RRFM_MOUSE	100.000	262	0	0	1	262	1	262	0.0	540
601886at33208_10090_0:000846	sp|Q9D8U7|DTWD1_MOUSE	100.000	304	0	0	1	304	1	304	0.0	636
412930at33208_10090_0:000915	sp|Q61712|DNJC1_MOUSE	100.000	552	0	0	1	552	1	552	0.0	1125
552949at33208_10090_0:000965	sp|Q9DAZ9|ANCHR_MOUSE	98.232	396	0	2	1	396	1	389	0.0	791
295556at33208_10090_0:0009e4	sp|Q8CC12|CDAN1_MOUSE	100.000	1239	0	0	1	1239	1	1239	0.0	2484
493505at33208_10090_0:000a10	sp|Q9WVQ5|MTNB_MOUSE	100.000	173	0	0	140	312	69	241	2.48e-128	365
553411at33208_10090_0:000ae9	sp|Q9DA19|CIR1_MOUSE	100.000	450	0	0	1	450	1	450	0.0	905
114954at33208_10090_0:000b1e	sp|P59114|PCIF1_MOUSE	100.000	706	0	0	1	706	1	706	0.0	1453
628231at33208_10090_0:000d21	sp|O88653|LTOR3_MOUSE	100.000	124	0	0	1	124	1	124	3.82e-88	251
625091at33208_10090_0:000ee4	sp|Q9CY16|RT28_MOUSE	100.000	186	0	0	1	186	1	186	3.06e-135	375
601875at33208_10090_0:000fad	sp|Q05CL8|LARP7_MOUSE	100.000	570	0	0	1	570	1	570	0.0	1157
601875at33208_10090_0:001497	sp|Q05CL8|LARP7_MOUSE	96.491	171	5	1	3	173	208	377	5.24e-104	308
203698at33208_10090_0:001633	sp|Q8VHK9|DHX36_MOUSE	28.643	796	442	31	237	992	189	898	1.69e-67	246
121127at33208_10090_0:00168e	sp|Q8BWZ3|NAA25_MOUSE	100.000	972	0	0	1	972	1	972	0.0	2012
584333at33208_10090_0:001b7f	sp|Q9D920|BORC5_MOUSE	100.000	176	0	0	21	196	20	195	1.80e-130	363
529708at33208_10090_0:001cdc	sp|Q8BI36|JKAMP_MOUSE	71.148	305	86	2	2	305	8	311	3.87e-160	448
613522at33208_10090_0:0020d1	sp|Q9EQI8|RM46_MOUSE	100.000	283	0	0	1	283	1	283	0.0	575
631767at33208_10090_0:00281b	sp|Q3UJ81|NEPR1_MOUSE	100.000	112	0	0	1	112	1	112	6.82e-82	235
524148at33208_10090_0:002a36	sp|Q8R3L2|TCF25_MOUSE	99.120	682	0	1	1	682	1	676	0.0	1390
534172at33208_10090_0:002cc0	sp|Q9DCU6|RM04_MOUSE	88.779	303	22	2	1	300	1	294	0.0	532
648300at33208_10090_0:0036e5	sp|Q9CR29|CCD43_MOUSE	100.000	222	0	0	1	222	1	222	2.52e-156	431
523376at33208_10090_0:0038d1	sp|A2A559|PGAP3_MOUSE	100.000	320	0	0	1	320	1	320	0.0	659
430701at33208_10090_0:003a81	sp|Q8QZR8|FA58B_MOUSE	100.000	250	0	0	1	250	1	250	0.0	520
590347at33208_10090_0:003d84	sp|Q8VDR5|TM267_MOUSE	100.000	215	0	0	24	238	1	215	5.30e-158	436
432287at33208_10090_0:003ded	sp|Q6NXW6|RAD17_MOUSE	100.000	688	0	0	1	688	1	688	0.0	1423
334272at33208_10090_0:003ff6	sp|Q7TN31|AGGF1_MOUSE	100.000	711	0	0	1	711	1	711	0.0	1467
641073at33208_10090_0:0042eb	sp|Q9QXX8|NUFP1_MOUSE	100.000	484	0	0	1	484	1	484	0.0	1006
393219at33208_10090_0:004611	sp|Q8C3I8|HGH1_MOUSE	100.000	393	0	0	1	393	1	393	0.0	778
508986at33208_10090_0:0047d6	sp|Q9JK23|PSMG1_MOUSE	100.000	289	0	0	1	289	1	289	0.0	595
206858at33208_10090_0:004813	sp|Q8K0F1|TBC23_MOUSE	97.854	699	0	1	1	699	1	684	0.0	1405
399439at33208_10090_0:004b63	sp|Q9EP82|WDR4_MOUSE	100.000	413	0	0	44	456	1	413	0.0	842
435253at33208_10090_0:004bb7	sp|Q6ZQL4|WDR43_MOUSE	100.000	677	0	0	1	677	1	677	0.0	1398
556461at33208_10090_0:004d69	sp|Q8VCR3|TM242_MOUSE	100.000	140	0	0	1	140	1	140	9.05e-99	279
569812at33208_10090_0:0052f6	sp|Q8BP78|F10C1_MOUSE	96.700	303	6	1	102	404	17	315	0.0	596
```

