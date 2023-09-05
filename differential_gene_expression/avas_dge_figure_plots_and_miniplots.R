# differential gene expression in A. vastus body vs developing tip
# main plot for figures 3-6
# by WRF 2022-12-04

library(dplyr)
library(ggplot2)

dge_data_file = "~/git/Aphrocallistes_vastus_genome/differential_gene_expression/Data/avas_v1_29_annotations_prot_fr_blastp_k10.annotation_table_clusters.tsv.gz"
dge_raw_counts_file = "~/git/Aphrocallistes_vastus_genome/differential_gene_expression/Data/Avas_DGE_quant.csv.gz"
dge_sample_info_file = "~/git/Aphrocallistes_vastus_genome/differential_gene_expression/Data/Avas_DGE_quant.info"

# read final DGE table
dge_data = read.table(dge_data_file, header=TRUE, sep="\t", quote = '"')
#head(dge_data)

dge_raw_counts = read.table(dge_raw_counts_file, header = TRUE, sep = "\t")
dge_sample_info = read.csv(dge_sample_info_file , stringsAsFactors = FALSE)
dge_sample_xpos = match( dge_sample_info$individum, sort(unique(dge_sample_info$individum) ) )

#sum(dge_data$log2FoldChange, na.rm = TRUE)

# easy test plot
# plot( dge_data$, -1*log10(dge_data$pvalue), 
#       xlim=c(-10,10),
#       pch=16, col=ifelse(is_p0001,"#045a8daa","#fe992922"))


# define arbitrary threshold of
PCUTOFF = 0.0001
is_sig_point = dge_data$padj < PCUTOFF

table( is_sig_point )
table( (dge_data$padj < PCUTOFF & dge_data$log2FoldChange >= 0) )
table( (dge_data$padj < PCUTOFF & dge_data$log2FoldChange < 0) ) 
dge_sig_file = "~/git/Aphrocallistes_vastus_genome/differential_gene_expression/Data/avas_v1_29_annotations_prot_fr_blastp_k10.annotation_table_clusters.p1e4.tsv"
dge_data_sig = filter(droplevels(dge_data), !is.na(dge_data$padj) & dge_data$padj < PCUTOFF)
write.table(dge_data_sig, dge_sig_file, sep="\t", quote=FALSE, row.names = FALSE)

is_sig_point_int = as.character(as.integer(is_sig_point) + 1)
is_sig_point_int[(dge_data$padj < PCUTOFF & dge_data$log2FoldChange < 0)] = "3"

fig3gg = ggplot(data = dge_data, aes( x=log2FoldChange, y=padj )) +
  theme(legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.title = element_text(size=23),
        plot.caption = element_text(hjust = 0, size=14) ) +
  labs(x="Log2-fold change", y="Adjusted p-value",
       caption="Supplemental Figure 09: Differential gene expression between body\nand osculum. With p-value threshold of 1e-4, 419 genes\nwere significantly differentially expressed,\n with 342 upregulated (purple) and 77 downregulated (blue).") +
  scale_x_continuous(limits=c(-10,10)) +
  scale_y_log10(limits=c(1e-101,1e5), breaks=c(1e-100, 1e-75, 1e-50, 1e-25), expand=c(0,0)) +
  scale_color_manual(values=c("#d8974f22", "#940c68aa", "#04688daa"), na.value = "#00000011") +
  geom_point( aes(color=is_sig_point_int), size=log10(dge_data$baseMean), shape=16) +
  #geom_point( color="#07cbe6bb", size=4, shape=16) +
  #annotate(geom="text", x=-Inf, y=Inf, size=5, hjust=0, vjust=1, label=paste(table(is_sig_point)[["TRUE"]], "significant out of", length(is_sig_point) ) ) +
  annotate(geom="text", x=0, y=Inf, size=4, hjust=0, vjust=1, label="  upregulated  >" ) +
  annotate(geom="text", x=0, y=Inf, size=4, hjust=1, vjust=1, label="<  downregulated  " ) +
  annotate(geom="segment", x=0, xend=0, y=10, yend=1e-100, size=1, alpha=0.3)
fig3gg
ggsave(filename = "~/git/Aphrocallistes_vastus_genome/differential_gene_expression/sfig9_padj.pdf", plot=fig3gg, device="pdf", width=7, height=7, paper="a4")


fig3gg = ggplot(data = dge_data, aes( x=log2FoldChange, y=baseMean )) +
  theme(legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.title = element_text(size=23),
        plot.caption = element_text(hjust = 0, size=14) ) +
  labs(x="Log2-fold change", y="Base mean counts",
       caption="Supplemental Figure 08: Differential gene expression between body\nand osculum. With p-value threshold of 1e-4, 419 genes\nwere significantly differentially expressed,\n with 342 upregulated (purple) and 77 downregulated (blue).") +
  scale_x_continuous(limits=c(-10,10)) +
  scale_y_log10(limits=c(0.1,1e5), expand=c(0,0)) +
  scale_color_manual(values=c("#d8974f22", "#940c68aa", "#04688daa"), na.value = "#00000011") +
  geom_point( aes(color=is_sig_point_int), size=2, shape=16) +
  annotate(geom="text", x=0, y=Inf, size=4, hjust=0, vjust=1, label="  upregulated  >" ) +
  annotate(geom="text", x=0, y=Inf, size=4, hjust=1, vjust=1, label="<  downregulated  " ) +
  annotate(geom="segment", x=0, xend=0, y=10, yend=1e-100, size=1, alpha=0.3)
fig3gg
ggsave(filename = "~/git/Aphrocallistes_vastus_genome/differential_gene_expression/sfig8_basemean.pdf", plot=fig3gg, device="pdf", width=7, height=7, paper="a4")



################################################################################


# mini plots for various loci figures
miniplot_target_names = c("Avas.s006.g409.i1", "Avas.s006.g479.i1", "Avas.s006.g1066.i1", "Avas.s006.g1068.i1", "Avas.s001.g969.i1", "Avas.s001.g970.i1",
                          "Avas.s015.g314.i1", "Avas.s003.g260.i1", # collagen hydroxylases
                          "Avas.s001.g887.i1", "Avas.s006.g483.i1", # proteases
                          "Avas.s014.g618.i1", "Avas.s014.g617.i1", # glassin
                          "Avas.s003.g337.i1", "Avas.s003.g335.i1", "Avas.s003.g338.i1", # TIL domains
                          "Avas.s006.g1096.i1", "Avas.s006.g1095.i1", "Avas.s015.g331.i1", "Avas.s015.g332.i1", # aquaporin
                          "Avas.s019.g225.i1", "Avas.s019.g227.i1", # DUOX
                          "Avas.s001.g1473.i1", "Avas.s001.g1475.i1", "Avas.s001.g1477.i1", "Avas.s001.g1480.i1", "Avas.s001.g1481.i1" # CA
)
miniplot_target_names = c( "Avas.s001.g887.i1", "Avas.s006.g483.i1", # proteases
                           "Avas.s014.g618.i1", # glassin
                           "Avas.s001.g1214.i1", "Avas.s014.g306.i1", "Avas.s009.g292.i1", # cathepsin-L-like
                           "Avas.s001.g39.i1", "Avas.s008.g653.i1", "Avas.s014.g146.i1", "Avas.s001.g983.i1", "Avas.s001.g47.i1", # cathepsin L
                           "Avas.s005.g1000.i1", "Avas.s005.g997.i1", "Avas.s004.g800.i1", "Avas.s019.g266.i1", "Avas.s019.g265.i1", "Avas.s001.g1617.i1", # cathepsin B
                           "Avas.s008.g921.i1", "Avas.s022.g202.i1", "Avas.s001.g1900.i1") # other cathepsins


miniplot_target_index = match( miniplot_target_names , dge_raw_counts$Name )

miniplot_x = rep(c(0.2, 0.6, 0.7, 0.8, 0.3, 0.4), each=4)
bt_order = sort(dge_sample_info$individum, index.return=TRUE)
pdf("~/git/Aphrocallistes_vastus_genome/differential_gene_expression/dge_gene_miniplots_v3.pdf", width=8, height=7)
par(mfrow = c(4,6), mar=c(3,3,3,1))
for (mti in miniplot_target_index){
  target_gene_raw = as.numeric(dge_raw_counts[mti,2:25])[bt_order$ix]
  bar_color = rep( c("#3075eedd", "#940c68dd"), each=12 )
  barplot(target_gene_raw, col=bar_color, main=gsub("Avas.", "", dge_raw_counts$Name[mti]), border = NA )
  text(0, 0.95*max(target_gene_raw), formatC(dge_data$padj[match(dge_raw_counts$Name[mti], dge_data$SeqName)], format="e", digits=1), pos=4)
  axis(1, at=c(7, 21), labels=c("Body", "Tip"), tick = FALSE, cex.axis=1)
}
dev.off()



################################################################################

# make Manhattan type plot of DGE across the large scaffolds

all2Ddata = read.table("~/git/Aphrocallistes_vastus_genome/synteny/vs_oopsacas/Avas.1.29_vs_oopsacas_gb.scaffold2d_points.tab", sep="\t", stringsAsFactors=FALSE)
categories = all2Ddata[,1]
is_scaf1 = which(categories=="s1")
scafdata1 = all2Ddata[is_scaf1,]
longscafs1 = c(0,scafdata1[1:50,6])
is_points = which(categories=="g")
pointsdata = droplevels(all2Ddata[is_points,])

#plot( all_avas_pointsdata[,6], -1*log10( (dge_data$padj)[points_dge_index]), pch=16, col="#d8974f22" )
#segments(longscafs1,0,longscafs1,100, lwd=0.5, col="#00000066")

has_avas_gene = !is.na(pointsdata[,3])
all_avas_pointsdata = pointsdata[has_avas_gene,]
points_dge_index = match(all_avas_pointsdata[,2], dge_data$SeqName )

avas_scaf_dge_df = data.frame(SeqName = all_avas_pointsdata[,2], 
                              position = all_avas_pointsdata[,6], 
                              baseMean = (dge_data$baseMean)[points_dge_index],
                              log2foldchange = (dge_data$log2FoldChange)[points_dge_index],
                              padj = -1*log10( (dge_data$padj)[points_dge_index]) )

manhplot_target_names = c("Avas.s006.g1066.i1", "Avas.s006.g1068.i1", "Avas.s001.g969.i1", "Avas.s001.g970.i1",
                          "Avas.s015.g314.i1", "Avas.s003.g260.i1", "Avas.s011.g507.i1", # collagen hydroxylases
                          "Avas.s001.g887.i1", "Avas.s006.g483.i1", # metallo proteases
                          "Avas.s014.g618.i1", # glassin
                          "Avas.s003.g337.i1", #"Avas.s003.g335.i1", "Avas.s003.g336.i1", # TIL domains
                          "Avas.s006.g1096.i1", # aquaporin
                          "Avas.s019.g225.i1", "Avas.s019.g227.i1", # DUOX and DUOXA
                          "Avas.s001.g868.i1", # bZIP TF
                          "Avas.s005.g703.i1", # the one actin out of 43
                          "Avas.s007.g465.i1", # cluster_00633, now hexaxilin
                          "Avas.s014.g646.i1", #protease
                          "Avas.s015.g491.i1" ) # cluster_00633, now hexaxilin

avas_scaf_dge_df_short = filter(avas_scaf_dge_df, SeqName %in% manhplot_target_names)
avas_scaf_dge_df_short
show_functions = c("bZIP", "protease", "COL", "COL", "P4H", "TIL", "actin", "protease", "COL", "COL", 
                   "aquaporin-9", "hexaxilin", "P3H", "glassin", "protease", "PLOD", "hexaxilin", "DUOX", "DUOXA")

point_color = rep("#00000011", length(avas_scaf_dge_df$padj))
is_pos_deg = avas_scaf_dge_df$padj > 4 & avas_scaf_dge_df$log2foldchange > 0
is_neg_deg = avas_scaf_dge_df$padj > 4 & avas_scaf_dge_df$log2foldchange < 0
point_color[is_pos_deg] = "#8d0468aa"
point_color[is_neg_deg] = "#04688daa"

ggman = ggplot(avas_scaf_dge_df, aes(x=position, y=padj)) +
  theme(legend.position="none",
        plot.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank() ) +
  labs(x=NULL, y="Adjusted p-value (-log)") +
  scale_x_continuous(breaks=(longscafs1[1:27] + (diff(longscafs1)/2)[1:27] ), labels=scafdata1[1:27,3], expand=c(0,0)) +
  scale_y_continuous(expand=c(0.01,0)) + 
  geom_point( color=point_color, size=log10(avas_scaf_dge_df$baseMean), shape=16) +
  annotate(geom="segment", x=longscafs1[1:28], y=rep(0,28), xend=longscafs1[1:28], yend=rep(100,28), alpha=0.3) +
  annotate(geom="text", x = avas_scaf_dge_df_short$position + 1e6, y = avas_scaf_dge_df_short$padj, label=show_functions, hjust=0)
ggman
ggsave(file="~/git/Aphrocallistes_vastus_genome/differential_gene_expression/dge_manhattan_plot.pdf", ggman, device = "pdf", width=8, height=3)
ggsave(file="~/git/Aphrocallistes_vastus_genome/differential_gene_expression/dge_manhattan_plot.png", ggman, device = "png", width=8, height=3, dpi=90)



#

