# trans-splicing processing steps #

Trans-splicing is detected in many animal phyla. Here, based on the mapping of the long reads to the genome, the "skip" `S` is extracted from the CIGAR string of each read in the SAM/BAM file using the script `get_read_skip_from_bam.py`. The skip is counted for all reads, and the leader sequence (the bases) is extracted for 42-45bp skips.

```
~/samtools-1.9/samtools view Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.sorted.bam | get_read_skip_from_bam.py - > Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.skip.tab
# Reading <stdin>, tracking S of 42-45bp
# Counted 25622489 lines for 25622489 reads
# 1152936 reads had either end matched to target leader length
# 57918 were likely polyA tails, and were ignored
```

The skipped sequences, putative leader sequences, are then counted. The forward ones end with `...AAATAC` and then would join `AG`. Reverse ones (i.e. those mapping to the reverse strand) are seen starting with `GTATTT...` or similar.

```
cut -f 4 Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.skip.tab | sort | uniq -c | sort -nr 
24531627 0
  51429 *
  11675 CTTAAAAAACTACTTTCAAAACTTCAAAACAAAACTAAATAC
   9192 GTATTTAGTTTTGTTTTGAAGTTTTGAAAGTAGTTTTTTAAG
   8574 GGGCTTTATAAAACTTCTACAAACTTCATACAAAACTAAATAC
   8526 GGGCTTATAAAACTACTTTACAAACTTCAAACAAAACTAAATAC
   8376 GGGCTTAAAAACTACTTTCAAAACTTCAAAACAAAACTAAATAC
   8315 GGGCTTTATAAAACTTCTACAAAACTTCATACAAATCTAAATAC
   7981 GGGGCTTAAAAACTACTTTCAAAACTTCAAAACAAAACTAAATAC
   7222 GGGCTTAAAAAACTACTTTCAAAACTTCAAACAAAACTAAATAC
   7115 GTATTTAGTTTTGTTTTGAAGTTTTGAAAGTAGTTTTTAAGCCC
   7056 GGGCTTTATAAAACTTCTACAAAACTTTATACAAATCTAAATAC
   7048 GTATTTAGTTTTGTTTTGAAGTTTTGAAAGTAGTTTTTAAGCCCC
   7020 GTATTTAGATTTGTATGAAGTTTTGTAGAAGTTTTATAAAGCCC
   6859 GTATTTAGTTTTGTATGAAGTTTGTAGAAGTTTTATAAAGCCC
   6849 GTATTTAGTTTTGTTTGAAGTTTGTAAAGTAGTTTTATAAGCCC
   6812 GGGCTTATAAAACTACTTTCAAAACTTCAAAACAAAACTAAATAC
   6484 CTTAAATAAACTACTTTACAAAACTTCAAACAAAACTAAATAC
   6085 GTATTTAGTTTTGTTTTGAAGTTTTGAAAGTAGTTTTATAAGCCC
   6073 GGGCTTTATAAAACTTCTACAAACTTCATACAAAAACTAAATAC
   5818 GGCTTAAAAAACTACTTTCAAAACTTCAAAACAAACTAAATAC
   5712 GCTTAAAAAACTACTTTCAAAACTTCAAAACAAAACTAAATAC

```

The sequences are aligned with MAFFT to determine the consensus.

`mafft leader_seq_scores.fasta > leader_seq_scores.aln`


