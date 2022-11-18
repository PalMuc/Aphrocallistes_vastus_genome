####################################
#
# Load required packages
#
#################################################################

library(reshape2)

####################################
#
# Functions
#
#################################################################

getDEGsofInterest <- function(cluster, DEGs) {
  
  DEGs_of_interest <- DEGs[DEGs$transcript %in% cluster,]
  DEGs_of_interest
  
}


make_barplot <- function(l, DEGs){
  
  list_to_use <- gene_lists[[l]]
  
  print(list_to_use)
  
  plot_name <- l
  
  if(length(list_to_use) > 1){
    
    DEGs <- all_DEGs
    
    x <- lapply(list_to_use, getDEGsofInterest, DEGs)
    
    melted_x <- melt(x, na.rm = T)
    df_x <- as.data.frame(t(acast(melted_x, variable ~ transcript)))
    df_x$transcript <- rownames(df_x)
    df_x <- merge(df_x, melted_x[,c(1,4)], by.y = "transcript")
    df_x <- df_x[!duplicated(df_x),]
    rownames(df_x) <- df_x$transcript
    df_x <- df_x[order(df_x$L1),]
    
    df_DEGs <- df_x[df_x$transcript %in% DEG_names$x,]
    
    df_DEGs$sign<-df_DEGs$log2FoldChange > 0
    
    #reorder levels to match graph to DEG table
    df_DEGs$transcript <- factor(df_DEGs$transcript, levels(as.factor(df_DEGs$transcript))[match(df_DEGs$transcript, levels(as.factor(df_DEGs$transcript)))])
    
    colorBars<-brewer.pal(5,"Set2")
    gLFC <- ggplot(df_DEGs) + geom_crossbar(aes(x=transcript, y=log2FoldChange, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE, fill=sign), fatten = 0.75, width=0.75, alpha=0.9) +
      scale_fill_manual(values=c(rgb(0,144,178,1,maxColorValue = 255),rgb(213,94,0,1,maxColorValue = 255)), name = "", labels = c("Downregulated","Upregulated") ) +
      scale_y_continuous(breaks = seq(-7,7,2), minor_breaks = seq(-2,6,2), limits = c(-7,7)) + 
      theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.grid.major = element_line(colour="black"), panel.grid.minor = element_line(colour="black", linetype = "dotted"), panel.grid.major.x = element_blank()) 
    
    
    ggsave(paste(plot_folder, plot_name, "_LFC_BarPlot.pdf", sep = ""), gLFC, device = "pdf", dpi = 300)


  }
}

###############################################
#
# out and plot folders
#
#########################################################################

plot_folder <- "../Plots/"


#################
#
# source external transcript lists -> this provides the necessary object gene_lists!!!
#
########################################

source("../Clusters/homolog_paralog_DEGs_makeLists.R")


#########################
#
# Read inputs
#
####################################################

DEG_names <- read.csv("../Outfiles/Avas_DEGs_FDR_0.05.names")

all_DEGs <- read.csv("../Outfiles/Avas_AllTranscripts_DESeq2Results.csv")
#DEGs <- read.csv("../Outfiles/Avas_DEGs_FDR_0.05.csv")
colnames(all_DEGs)[1] <- "transcript"


#####################################

lapply(names(gene_lists), make_barplot, all_DEGs)



