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
getDEGsofInterest <- function(cluster, DEG_df) {
  
  DEGs_of_interest <- DEG_df[DEG_df$transcript %in% cluster,]
  DEGs_of_interest
}


make_box_plot <- function(l, DEGs){
  
  list_to_use <- gene_lists[[l]]
  
  if(length(list_to_use) > 1){
    
    plot_name <- l
    
    filtered_orthologous_clusters <- list_to_use[lapply(list_to_use, length) > 1]
    
    x <- lapply(filtered_orthologous_clusters, getDEGsofInterest, DEGs)
    
    if(length(x) > 0){
      melted_x <- melt(x, na.rm = T)
      
      df_x <- as.data.frame(t(acast(melted_x, variable ~ transcript)))
      
      df_x$transcript <- rownames(df_x)
      
      df_x <- merge(df_x, melted_x[,c(1,4)], by.y = "transcript")
      
      df_x <- df_x[!duplicated(df_x),]
      
      g<-ggplot(df_x, aes(x=L1, y=log2FoldChange)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(size=(log2FoldChange)^2, color=log2FoldChange), width = 0.1) +
        scale_colour_gradient2(low = "darkblue", high = "darkorange", mid = "lightgrey", midpoint = 0, na.value = NA) +
        theme_bw() + 
        #theme(axis.text.x = element_blank(),panel.grid.major = element_line(colour="black"), panel.grid.minor = element_line(colour="black", linetype = "dotted"), panel.grid.major.x = element_blank())  
        theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.grid.major = element_line(colour="black"), panel.grid.minor = element_line(colour="black", linetype = "dotted"), panel.grid.major.x = element_blank()) 
      
      ggsave(paste(plot_folder, plot_name, "_LFC_BoxPlot.pdf", sep = ""), g, device = "pdf", dpi = 300)
      
    }
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
all_DEGs <- read.csv("../Outfiles/Avas_AllTranscripts_DESeq2Results.csv")
#DEGs <- read.csv("../Outfiles/Avas_DEGs_FDR_0.05.csv")
colnames(all_DEGs)[1] <- "transcript"


#####################################
#
# gen box plots
#
###############################################################

lapply(names(gene_lists), make_box_plot, all_DEGs)


