####################################
#
# Load required packages
#
#################################################################


library(pheatmap)
library(reshape2)


####################################
#
# Functions
#
#################################################################
zscore_transform <- function(count_df){
  
  row_mean <- apply(count_df,1,mean)
  row_sd <- apply(count_df,1,sd)
  
  (count_df-row_mean)/row_sd
  
}


getCountssofInterest <- function(cluster, count_matrix) {
  
  count_matrix_of_interest <- count_matrix[rownames(count_matrix) %in% cluster,]
  count_matrix_of_interest

}

make_heat_plot <- function(l, row_labels, t_counts){
  
  list_to_use <- gene_lists[[l]]
  
  plot_name <- l
  
  x <- lapply(list_to_use, getCountssofInterest, t_counts)
  
  melted_x <- melt(x, na.rm = T)
  df_x <- as.data.frame(t(acast(melted_x, variable ~ transcript)))
  df_x$transcript <- rownames(df_x)
  df_x <- merge(df_x, melted_x[,c(1,4)], by.y = "transcript")
  df_x <- df_x[!duplicated(df_x),]
  rownames(df_x) <- df_x$transcript
  df_x <- df_x[order(df_x$L1),]
  
  df_DEGs <- df_x[df_x$transcript %in% DEG_names$x,]
  
  if(nrow(df_DEGs) > 2) {
    col_annotation_df <- data.frame(Body_Part = factor(rep(c("Body","Tip"), c(3,3)), ordered = T), row.names = colnames(df_x[,2:7]))
    row_annotation_df <- data.frame(Cluster = factor(df_DEGs[,8], ordered = T), row.names = df_DEGs[,1])
    pheatmap(df_DEGs[,2:7], show_rownames = row_labels, cluster_rows = T, annotation_col = col_annotation_df, cluster_cols = T, main=plot_name, filename = paste(plot_folder, plot_name,"_heatplot",".pdf",sep = ""))
  }
  
  write.csv(df_DEGs[,c(1,8)],paste(outfile_folder,plot_name,"_DEGs.csv", sep = ""),row.names = F)
  
} 

###############################################
#
# out and plot folders
#
#########################################################################

plot_folder <- "../Plots/"
outfile_folder <- "../Outfiles/"

#################
#
# source external transcript lists -> this provides the necessary object gene_lists!!!
#
########################################

source("../Clusters/homolog_paralog_DEGs_makeLists.R")

#####################################

#########################
#
# Read inputs
#
####################################################

DEG_names <- read.csv("../Outfiles/Avas_DEGs_FDR_0.05.names")
avas_counts <- read.csv("../Outfiles/Avas_DESeq2_counts.csv")

# Z-score transform and plot heatplot

avas_counts[,2:7] <- zscore_transform(avas_counts[,2:7])
rownames(avas_counts) <- avas_counts[,1]
colnames(avas_counts)[1] <- "transcript"

#FALSE = no transcript names in plot
#TRUE = transcript names in plot
lapply(names(gene_lists), make_heat_plot, FALSE, avas_counts)

