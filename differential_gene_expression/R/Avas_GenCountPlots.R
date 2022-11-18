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


make_box_plot <- function(l, counts, DEGs){
  
  list_to_use <- gene_lists[[l]]
  plot_name <- l
  print(l)
  x <- lapply(list_to_use, getDEGsofInterest, counts)
    
  if(length(x) > 0){
    melted_x <- melt(x, na.rm = T)
    
    melted_x[grep("B", melted_x$variable),5] <- "Body"
    melted_x[grep("B", melted_x$variable, invert = T),5] <- "Tip"
      
    deg_melted_x<-melted_x[melted_x$transcript %in% DEGs$transcript,]
      
    g<-ggplot(deg_melted_x, aes(x=V5, y=value)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(size=(value)^2, color=value), width = 0.1) +
      facet_wrap(vars(transcript), scales = "free_y")  +
      scale_colour_gradient2(low = "darkblue", high = "darkorange", mid = "lightgrey", midpoint = 0, na.value = NA) +
      theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.grid.major = element_line(colour="black"), panel.grid.minor = element_line(colour="black", linetype = "dotted"), panel.grid.major.x = element_blank()) 
      
    ggsave(paste(plot_folder, plot_name, "_Count_BoxPlot.pdf", sep = ""), g, device = "pdf", dpi = 300)
      
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
transcript_counts <- read.csv("../Outfiles/Avas_DESeq2_counts.csv")

DEG_names <- read.csv("../Outfiles/Avas_DEGs_FDR_0.05.names")
colnames(transcript_counts)[1] <- "transcript"
colnames(DEG_names)[1] <- "transcript"

#####################################
#
# gen box plots
#
###############################################################

lapply(names(gene_lists), make_box_plot, transcript_counts, DEG_names)


