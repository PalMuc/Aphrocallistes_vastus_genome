# plot histogram of read skips from long RNAseq mapping
# created by WRF  last modified 2022-12-14

# ~/samtools-1.9/samtools view Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.sorted.bam | get_read_skip_from_bam.py - > Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.skip.tab

skipdata_file = "~/genomes/aphrocallistes_vastus_PORI/Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.skip.tab.gz"
skipdata = read.table(skipdata_file, sep="\t")

forward_skips = table(skipdata[,2])
reverse_skips = table(skipdata[,3])

pdf(file="~/genomes/aphrocallistes_vastus_PORI/Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.skip.pdf", height=6, width=7)
plot(0,0,type='n', xlim=c(0,100), ylim=c(0,500000), ylab="Count", xlab="S length (bp)", main=basename(skipdata_file) )
lines( as.integer(unlist(dimnames(forward_skips))), forward_skips )
lines( as.integer(unlist(dimnames(reverse_skips))), reverse_skips )
text( as.integer(unlist(dimnames(forward_skips))), forward_skips, as.integer(unlist(dimnames(forward_skips))), pos=4)
dev.off()






#