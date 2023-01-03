# replot single cell data
# last modified 2022-12-23

library(dplyr)
library(ggplot2)

################################################################################
# analyze Sebe-Pedros 2018 Amphimedon data

aque_exp_data_file = "~/project/supplements/sebe-pedros_2018/aq_adult_gene_sums.tab"
#aque_exp_data_file = "~/project/supplements/sebe-pedros_2018/aq_larva_gene_sums.tab"
aque_exp_data = read.table(aque_exp_data_file,header=FALSE, sep="\t")
head(aque_exp_data)
dim(aque_exp_data)

aque_all_genes_expr = aque_exp_data$V2
is_zero_expr = (aque_all_genes_expr==0)

aque_cathepsin_names = c("Aqu2.1.41049_001", "Aqu2.1.43732_001", "Aqu2.1.43733_001", "Aqu2.1.22535_001", "Aqu2.1.39162_001")
aque_cathepsin_index = match(aque_cathepsin_names,aque_exp_data$V1)
aque_silicatein_names = c("Aqu2.1.42494_001", "Aqu2.1.42495_001", "Aqu2.1.41046_001", "Aqu2.1.41047_001", "Aqu2.1.41048_001", "Aqu2.1.29399_001")
aque_silicatein_index = match(aque_silicatein_names,aque_exp_data$V1)



pdf(file="~/git/Aphrocallistes_vastus_genome/biomineralization/aq_adult_gene_sums.pdf", width=6, height=6, paper="a4")
#pdf(file="~/git/Aphrocallistes_vastus_genome/biomineralization/aq_larva_gene_sums.pdf", width=6, height=6, paper="a4")
par(mar=c(10,4.5,2,2))
plot( 1:length(aque_all_genes_expr), aque_all_genes_expr, type='n', frame.plot=FALSE,
      xlab="Ordinal rank of 49905 A. queenslandica genes", ylab="Gene expression (sum of 4966 cells)",
      ylim=c(0,50000) )
lines( (1:length(aque_all_genes_expr))[!is_zero_expr], aque_all_genes_expr[!is_zero_expr], col="#08306b", lwd=3)
lines( (1:length(aque_all_genes_expr))[is_zero_expr], aque_all_genes_expr[is_zero_expr], col="#888888", lwd=3)
points(aque_silicatein_index, aque_all_genes_expr[aque_silicatein_index] - 1000, pch=17, col="#9ecae1", lwd=2)
points(aque_cathepsin_index, aque_all_genes_expr[aque_cathepsin_index] - 1000, pch=17, col="#fe9929", lwd=2)
segments( max(which(aque_all_genes_expr>4966)), 0, y1=300000, col="#00000088", lty=2 )
segments( 0, 4966, x1=50000, col="#00000088", lty=2 )
text(650,6000,"Expr > 5000 total", pos=4)
legend(20000, 15000, legend=c("Silicateins","Cathepsins"), pch=17, col=c("#9ecae1", "#fe9929"), bty="n")
mtext("Supplemental Figure 11:\nPlot of A. queenslandica body single cell RNAseq\nfrom Sebe-Pedros et al (2018).", 1, line=8, cex = 1.3, font=2)
dev.off()


################################################################################
# analyze Musser 2021 Spongilla data

spongilla_data_file = "~/project/musser2021_sponge-single-cell/supplementary_files/Data_S1_sheet3_diff_exp_cell_types.tab"
spongilla_data = read.table(spongilla_data_file, header=TRUE, sep="\t", quote='"', stringsAsFactors = FALSE)
names(spongilla_data)
table(spongilla_data$Cell.Type)
# Amb apnPin1 apnPin2     Apo     Arc  basPin    Chb1    Chb2     Cho     Grl incPin1 incPin2     Lph 
# 921    2140    1374    3062    6153    1248    3356    2183    2501     699    1078    1175     793 
# Mes1    Mes2    Mes3    Met1    Met2    Myp1    Myp2     Nrd     Scl     Scp 
# 400     356     485    1603    1856    3308    2294     780    1654    1073 

sclerocyte_data = filter(spongilla_data, Cell.Type=="Scl")
#sclerocyte_genes = filter(sclerocyte_data, Mean.Log.Fold.Diff. >= 1)
sclerocyte_genes = filter(sclerocyte_data, Percent.cells.expressing..in.cluster.-Percent.cells.expressing..outside.cluster. >= 0.4)

ggplot(sclerocyte_data, aes(x=Percent.cells.expressing..outside.cluster., y=Mean.Log.Fold.Diff.)) +
  theme(legend.position="none") +
  scale_x_continuous(limits=c(0,0.5)) +
  scale_color_continuous(type="gradient", high="#dddddd", low="#67001f") +
  geom_point(aes(color=Percent.cells.expressing..in.cluster.), size=3, alpha=0.8) +
  annotate(geom="text", x=sclerocyte_genes$Percent.cells.expressing..outside.cluster., 
           y=sclerocyte_genes$Mean.Log.Fold.Diff. , 
           label=sapply(strsplit( as.character(sclerocyte_genes$Automated.Gene.Name..in.seurat.object.and.some.suppl..Figs..) ," "), `[`, 1),
           hjust=0)

ggplot(sclerocyte_data, aes(x=Percent.cells.expressing..outside.cluster., y=Percent.cells.expressing..in.cluster.)) +
  theme(legend.position="none") +
  scale_x_continuous(limits=c(0,1), expand = c(0.02,0) ) +
  scale_y_continuous(limits=c(0,1), expand = c(0.02,0) ) +
  scale_color_continuous(type="gradient", high="#dddddd", low="#67001f") +
  geom_point(aes(size=Mean.Log.Fold.Diff.), color="#67001f", alpha=0.5) +
  annotate(geom="text", x=sclerocyte_genes$Percent.cells.expressing..outside.cluster., 
           y=sclerocyte_genes$Percent.cells.expressing..in.cluster., 
           label=sapply(strsplit( as.character(sclerocyte_genes$Automated.Gene.Name..in.seurat.object.and.some.suppl..Figs..) ," "), `[`, 1),
           hjust=0)


#