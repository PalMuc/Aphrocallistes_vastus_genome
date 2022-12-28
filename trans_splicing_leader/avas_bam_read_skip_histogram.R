# plot histogram of read skips from long RNAseq mapping
# created by WRF  last modified 2022-12-20

# generated with command:
# ~/samtools-1.9/samtools view Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.sorted.bam | get_read_skip_from_bam.py - > Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.skip.tab

skipdata_file = "~/git/Aphrocallistes_vastus_genome/trans_splicing_leader/Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.skip.tab.gz"
skipdata = read.table(skipdata_file, sep="\t")

forward_skips = table(skipdata[,2])
reverse_skips = table(skipdata[,3])

pdf(file="~/git/Aphrocallistes_vastus_genome/trans_splicing_leader/Aphrocallistes_instagraal_l4_n50_c1_N5.polished_sorted_hardmasked.ONT-RNA_minimap2.skip.pdf", 
    height=7, width=7, paper="a4")
par(mar=c(10,4.5,1,1))
plot(0,0,type='n', xlim=c(0,100), ylim=c(0,500000), 
     ylab="Count", xlab="S length (bp)", 
     frame.plot=FALSE)
     #main=basename(skipdata_file) )
text( as.integer(unlist(dimnames(forward_skips))), forward_skips, as.integer(unlist(dimnames(forward_skips))), pos=4)
lines( as.integer(unlist(dimnames(forward_skips))), forward_skips, lwd=2, col="#0101a9dd" )
lines( as.integer(unlist(dimnames(reverse_skips))), reverse_skips, lwd=2, col="#a4a6dedd" )
legend("topright", legend=c("forward","reverse"), lwd=5, col=c("#0101a9dd","#a4a6dedd"), bty="n")
mtext("Supplemental Figure 4:\nHistogram of read skip lengths from long RNA reads.", 1, line=8, cex = 1.4, font=2)
dev.off()






#