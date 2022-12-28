# draw synteny schematic for Avastus vs Ominuta
# created by WRF 2022-10-24
# last modified 2022-12-12

setwd("~/git/Aphrocallistes_vastus_genome/")

all2Ddata = read.table("synteny/vs_oopsacas/Avas.1.29_vs_oopsacas_gb.scaffold2d_points.tab.gz", sep="\t", stringsAsFactors=FALSE)
categories = all2Ddata[,1]
is_scaf1 = which(categories=="s1")
scafdata1 = all2Ddata[is_scaf1,]
longestscaf1 = max(scafdata1[,6])
is_longscafs1 = which(as.numeric(scafdata1[,5]) > 0.01)
longscafs1 = c(0, scafdata1[,6][is_longscafs1] )
is_scaf2 = which(categories=="s2")
scafdata2 = all2Ddata[is_scaf2,]
longestscaf2 = max(scafdata2[,6])
is_longscafs2 =  which(as.numeric(scafdata2[,5]) > 0.0009)
longscafs2 = c(0, scafdata2[,6][is_longscafs2] )
is_points = which(categories=="g")
pointsdata = droplevels(all2Ddata[is_points,])
longscafpoints = pointsdata[,6] < 270000000
genome_x = pointsdata[,7][longscafpoints]
genome_y = pointsdata[,6][longscafpoints]
xmax = tail( pretty(longscafs2), n=1)
ymax = tail( pretty(longscafs1), n=1)
longscafs_x = longscafs2
longscafs_y = longscafs1
nscafs_x = length(longscafs2)
nscafs_y = length(longscafs1)
xmax_mb = round(xmax / 1000000)
ymax_mb = round(ymax / 1000000)
pointsize = log10(as.numeric(pointsdata[,8])) / 6

# set scaffolds, which will be reset again for the lower panel
allg1chr_genes = droplevels(pointsdata[pointsdata[,3]=="Aphrocallistes_vastus_HiC-scaffold_022" & !is.na(pointsdata[,3]),])
table(allg1chr_genes[,5],useNA = "ifany")
sum(table(allg1chr_genes[,5],useNA = "ifany"))
allg2chr_genes = pointsdata[pointsdata[,5]=="JAKMXF010000133.1" & !is.na(pointsdata[,5]),]
selectedpoints = allg1chr_genes[allg1chr_genes[,5]=="JAKMXF010000133.1",]
is_syntentic1_block = ( allg1chr_genes[,5]=="JAKMXF010000133.1" )
table(is_syntentic1_block, useNA = "ifany")
is_syntentic2_block = ( allg2chr_genes[,3]=="Aphrocallistes_vastus_HiC-scaffold_022" )
unique(allg2chr_genes[is_syntentic2_block,4])
table(is_syntentic2_block, useNA = "ifany")
synteny1_color = ifelse(is_syntentic1_block, "#0e2f9988", "#d76404cc")
synteny1_color[is.na(synteny1_color)] = "#00000088"
synteny2_color = ifelse(is_syntentic2_block, "#0e2f9988", "#d76404cc")
synteny2_color[is.na(synteny2_color)] = "#00000088"

# get names of genes in the cluster, to color separately
cluster1_names = read.table("synteny/vs_oopsacas/avas_s022_tandem_cluster_1.names")
cluster2_names = read.table("synteny/vs_oopsacas/avas_s022_tandem_cluster_2.names")
#cluster3_names = read.table("synteny/vs_oopsacas/avas_s022_tandem_cluster_3.names")
# table(allg1chr_genes[,2] %in% cluster1_names[,1])
# table(allg1chr_genes[,2] %in% cluster2_names[,1])
# table(allg1chr_genes[,2] %in% cluster3_names[,1])
synteny1_color[allg1chr_genes[,2] %in% cluster1_names[,1]] = "#41ca3ccc"
synteny1_color[allg1chr_genes[,2] %in% cluster2_names[,1]] = "#ae3786cc"
#synteny1_color[allg1chr_genes[,2] %in% cluster3_names[,1]] = "#37720ccc"
synteny2_color[allg2chr_genes[,2] %in% gsub("O_minuta\\|","",cluster1_names[,1])] = "#41ca3ccc"
synteny2_color[allg2chr_genes[,2] %in% gsub("O_minuta\\|","",cluster2_names[,1])] = "#ae3786cc"
#synteny2_color[allg2chr_genes[,2] %in% gsub("O_minuta\\|","",cluster3_names[,1])] = "#37720ccc"
 

#
outputfile = "synteny/vs_oopsacas/Avas.v1.29_vs_oopsacas_gb.scaf_22_synteny.pdf"
#outputfile = "synteny/vs_oopsacas/Avas.v1.29_vs_oopsacas_gb.scaf_21_22_synteny.pdf"
pdf(file=outputfile, width=8, height=4)
#par( mfrow=c(2,1) )
par( mar=c(3,1,2,1) )
# empty plot
plot(0,0,type="n", xlim=c(-500,1800), ylim=c(-1,1), xlab="", ylab="", 
     frame.plot=FALSE, axes=FALSE)

# draw Avas scaf 22
rect(0, 0.4, (longscafs1[23]-longscafs1[22])/1000, 0.6, col="#c0c0cc")
segments( (allg1chr_genes[,6]-longscafs1[22])/1000, rep(0.38,dim(allg1chr_genes)[1]), 
          (allg1chr_genes[,6]-longscafs1[22])/1000, rep(0.62,dim(allg1chr_genes)[1]), 
          lwd=2, col=synteny1_color )
text(-500, 0.5, "A. vastus \nscaffold_022", cex=1.1, pos=4)

# draw Omin scaf 133
max_scaf_length = (longscafs2[22]-longscafs2[21])/1000
center_offset = 465
rect(center_offset, -0.4, center_offset + (longscafs2[22]-longscafs2[21])/1000, -0.6, 
     col="#c0c0cc")
segments( center_offset + max_scaf_length-(allg2chr_genes[,7]-longscafs2[21])/1000, rep(-0.62,dim(allg2chr_genes)[1]), 
          center_offset + max_scaf_length-(allg2chr_genes[,7]-longscafs2[21])/1000, rep(-0.38,dim(allg2chr_genes)[1]), 
          lwd=2, col=synteny2_color )
text(-500, -0.5, "O. minuta \nJAKMXF010000133.1 \n(reversed)", cex=1.1, pos=4)

linker_colors = rep("#00006655",nrow(selectedpoints))
linker_colors[selectedpoints[,2] %in% cluster1_names[,1]] = "#41ca3c66"
linker_colors[selectedpoints[,2] %in% cluster2_names[,1]] = "#ae378666"
#linker_colors[selectedpoints[,2] %in% cluster3_names[,1]] = "#37720c66"

segments( (selectedpoints[,6]-longscafs1[22])/1000, rep(0.38,dim(selectedpoints)[1]), 
          center_offset + max_scaf_length-(selectedpoints[,7]-longscafs2[21])/1000, rep(-0.38,dim(selectedpoints)[1]), 
          col=linker_colors, lwd=1.0)
mtext("Length (Mb)",side=1, cex=1.3, line=1, at = 900)
axis(3,at=seq(0,1800,200), labels = seq(0,1.8,0.2), cex.axis=1.2,line=-2)
axis(1,at=center_offset + seq(0,800,200), labels = seq(0,0.8,0.2), cex.axis=1.2,line=-2)
dev.off()



# set scaffolds, which were previously set for the upper panel
allg1chr_genes = pointsdata[pointsdata[,3]=="Aphrocallistes_vastus_HiC-scaffold_021",]
allg2chr_genes = pointsdata[pointsdata[,5]=="JAKMXF010000111.1",]
selectedpoints = allg1chr_genes[allg1chr_genes[,5]=="JAKMXF010000111.1",]
is_syntentic1_block = ( allg1chr_genes[,5]=="JAKMXF010000111.1" )
is_syntentic2_block = ( allg2chr_genes[,3]=="Aphrocallistes_vastus_HiC-scaffold_021" )
synteny1_color = ifelse(is_syntentic1_block, "#0e2f9988", "#d76404cc")
synteny1_color[is.na(synteny1_color)] = "#00000088"
synteny2_color = ifelse(is_syntentic2_block, "#0e2f9988", "#d76404cc")
synteny2_color[is.na(synteny2_color)] = "#00000088"


outputfile = "synteny/vs_oopsacas/Avas.v1.29_vs_oopsacas_gb.scaf_21_synteny.pdf"
pdf(file=outputfile, width=8, height=5)
par( mar=c(3,1,2,1) )

plot(0,0,type="n", xlim=c(-500,1800), ylim=c(-1,1), xlab="", ylab="", 
     frame.plot=FALSE, axes=FALSE)

# draw Avas scaf 21
rect(0, 0.4, (longscafs1[22]-longscafs1[21])/1000, 0.6, col="#c0c0cc")
segments( (allg1chr_genes[,6]-longscafs1[21])/1000, rep(0.38,dim(allg1chr_genes)[1]), 
          (allg1chr_genes[,6]-longscafs1[21])/1000, rep(0.62,dim(allg1chr_genes)[1]), 
          lwd=2, col=synteny1_color )
text(-500, 0.5, "A. vastus \nscaffold_021", cex=1.1, pos=4)

# draw Omin scaf 111
center_offset = 0
rect(center_offset, -0.4, center_offset + (longscafs2[2]-longscafs2[1])/1000, -0.6, 
     col="#c0c0cc")
segments( center_offset + (allg2chr_genes[,7]-longscafs2[1])/1000, rep(-0.62,dim(allg2chr_genes)[1]), 
          center_offset + (allg2chr_genes[,7]-longscafs2[1])/1000, rep(-0.38,dim(allg2chr_genes)[1]), 
          lwd=2, col=synteny2_color )
text(-500, -0.35, "O. minuta \nJAKMXF010000111.1 \n(reversed)", cex=1.1, pos=4)
segments( (selectedpoints[,6]-longscafs1[21])/1000, rep(0.38,dim(selectedpoints)[1]), 
          center_offset + (selectedpoints[,7]-longscafs2[1])/1000, rep(-0.38,dim(selectedpoints)[1]), 
          col="#00000066", lwd=1.0)
mtext("Length (Mb)",side=1, cex=1.3, line=1, at = 900)
axis(3,at=seq(0,1800,200), labels = seq(0,1.8,0.2), cex.axis=1.2,line=-2)
axis(1,at=seq(0,1800,200), labels = seq(0,1.8,0.2), cex.axis=1.2,line=-2)
dev.off()





################################################################################
# FIGURE 2 - exons per gene on scaffold 022

av_sc22_exons = read.table("genomic_overview/Avas.v1.29_annotations.fr.scaf_022.exon_hist.tab")
av_sc22_exons.sorted = sort(av_sc22_exons[,1],index.return=TRUE)
om_sc133_exons = read.table("genomic_overview/JAKMXF01.1.gbff.scaf_133.exon_hist.tab")
om_sc133_exons.sorted = sort(om_sc133_exons[,1],index.return=TRUE)
pdf("Avas_scaf_022_Omin_scaf_133.exon_hist.pdf", width = 4, height = 3, useDingbats = FALSE)
par(mar=c(4,4,1,0))
plot( av_sc22_exons.sorted$x - 0.2, av_sc22_exons[av_sc22_exons.sorted$ix,2], type='h',
      xlab ="", ylab = "", frame.plot = FALSE , lwd=3, col="#0e2f99cc",
      cex.axis=1.2, cex.lab=1.4, xlim = c(0,70))
segments( om_sc133_exons.sorted$x+0.2, rep(0,nrow(om_sc133_exons)), om_sc133_exons.sorted$x+0.2, 
          om_sc133_exons[ om_sc133_exons.sorted$ix,2], lwd = 3,
          col =  "#d76404cc")
legend(5,300,legend = c("A. vastus", "O. minuta"), pch = 15, col = c("#0e2f99cc", "#d76404cc"), bty = "n",
       text.font = 3, pt.cex = 2)
legend(25,300,legend = c("- scaffold 022", "- scaffold 133"), bty = "n")
text(58,50,"USP24")
segments(52,10,54,30, col = "#000000bb", lwd = 0.5)
segments(63,10,61,30, col = "#000000bb", lwd = 0.5)
mtext("Exons per gene", 1, line = 2.4, cex = 1.4)
mtext("N genes", 2, line = 2.4, cex = 1.4)
dev.off()






################################################################################
# FIGURE 2 - intron lengths on scaffold 022

avas_intron_file = "~/git/Aphrocallistes_vastus_genome/genomic_overview/Avas.v1.29_annotations.fr.introns.txt.gz"
avas_intron = read.table(avas_intron_file, sep = "\t")

omin_intron_file = "~/git/Aphrocallistes_vastus_genome/genomic_overview/JAKMXF01.1.gbff.introns.txt.gz"
omin_intron = read.table(omin_intron_file, sep = "\t")

h1 = hist(avas_intron$V5, breaks=c(seq(0,5000,50),1000000), plot = FALSE )
h2 = hist(omin_intron$V5, breaks=c(seq(0,5000,50),1000000), plot = FALSE )
plot(0,0, type='n', xlim=c(0,1500), ylim=c(0,10000),
     xlab="Intron length (bp)", ylab="Introns (n)")
lines( h1$mids, h1$counts, lwd=4, col="#0e2f9988")
lines( h2$mids, h2$counts, lwd=4, col="#d7640488")

avas_s022_introns = avas_intron$V5[avas_intron$V1=="Aphrocallistes_vastus_HiC-scaffold_022"]
omin_s133_introns = omin_intron$V5[omin_intron$V1=="JAKMXF010000133.1"]
h1 = hist(avas_s022_introns, breaks=c(seq(0,5000,50),1000000), plot = FALSE )
h2 = hist(omin_s133_introns, breaks=c(seq(0,5000,50),1000000), plot = FALSE )
plot(0,0, type='n', xlim=c(0,1500), ylim=c(0,360),
     xlab="Intron length (bp)", ylab="Introns (n)")
segments( h1$mids-5, rep(0,length(h1$counts)), h1$mids-5, h1$counts, lwd=4, col="#0e2f9988")
segments( h2$mids+5, rep(0,length(h2$counts)), h2$mids+5, h2$counts, lwd=4, col="#d7640488")

# significance test, since reviewers senselessly ask for this
combined_target_scaf_introns = data.frame(hists=c(h1$counts,h2$counts),
                                          sp=c(rep("Av",length(h1$counts)),rep("Om",length(h2$counts))) )
wilcox.test(formula=hists ~ sp, data=combined_target_scaf_introns)
# Wilcoxon rank sum test with continuity correction
# 
# data:  hists by sp
# W = 6542.5, p-value = 0.0001489
# alternative hypothesis: true location shift is not equal to 0

avas_s022_introns = avas_intron$V5[avas_intron$V1=="Aphrocallistes_vastus_HiC-scaffold_022"]
combined_target_scaf_introns = data.frame(hists=c(avas_s022_introns, omin_s133_introns),
                                          sp=c( rep("Av",length(avas_s022_introns))
                                                ,rep("Om",length(omin_s133_introns)) ) )
wilcox.test(formula=hists ~ sp, data=combined_target_scaf_introns)


s022_bp_sum = sum(avas_s022_introns)
s133_bp_sum = sum(omin_s133_introns)
s022_bp_diff = round( (s022_bp_sum - s133_bp_sum)/1000 )

pdf(file = "avas_s022_introns_vs_omin_s133_v1.pdf", width = 4, height = 4, useDingbats = FALSE)
par(mar=c(4.1,4.1,1,1))
plot(0,0, type='n', xlim=c(0,810), ylim=c(0,5000),
     xlab="", ylab="",
     cex.axis=1.3, cex.lab=1.3, frame.plot=FALSE)
mtext("Ordinal rank (long to short)", 1, line = 2.4, cex = 1.4)
mtext("Intron length (bp)", 2, line = 2.4, cex = 1.4)
polygon( c( seq(1,length(avas_s022_introns),1),length(avas_s022_introns),0 ), 
         c(sort(avas_s022_introns,decreasing = TRUE),0,0), col="#0e2f99cc", border = NA)
polygon( c( seq(1,length(omin_s133_introns),1),length(omin_s133_introns),0 ), 
         c(sort(omin_s133_introns,decreasing = TRUE),0,0), col="#ec8732ff", border = NA)
legend(20,5000,legend = c("A. vastus", "O. minuta"), pch = 15, col = c("#0e2f99cc", "#d76404cc"), bty = "n",
       text.font = 3, pt.cex = 2)
legend(255,5000,bty = "n", legend = c("- scaffold 022","- scaffold 133") ) 
#legend(250,5000,bty = "n", legend = c(paste("scaffold 022 :", round( s022_bp_sum/1000 ),"kbp"),
#                  paste("scaffold 133 :", round( s133_bp_sum/1000 ),"kbp") ) )
text(110,1200,paste(s022_bp_diff,"intron kbp difference"), pos=4 )
segments(90,600,120,1000, col = "#000000bb", lwd = 0.7)
dev.off()



#