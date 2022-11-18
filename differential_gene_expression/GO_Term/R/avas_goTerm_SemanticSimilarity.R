###
#
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("topGO")
#
#BiocManager::install("GOSemSim")

install.packages("../org.Avastus.eg.db/", repos = NULL)
library(topGO)
library(GOSemSim)
library(dendextend)
library(circlize)
library(RColorBrewer)
library(igraph)
library(grid)
library(gridExtra)
library(gtable)
library(dynamicTreeCut)
library(WGCNA)

setwd("/netvolumes/srva229/home/svargas/Checked_Clean_Stuff/Aphrocallistes_vastus_genome/05-ANALYSES/differential_gene_expression/GO_Term/R/")

####################################
#
# set files paths
#
###########################################################

#change this to use MF or CC ontologies and modify the TopGO section below accordingly
#once the wd is defined you can point to the files using paths relative to the wd
#Biological process table
processPath <- "../GOTerms/avas.BP.gos.4R.clean.topGO"

#deg list
deg_namesPath <- "../GOIs/Avas_downDEGs_FDR_0.05_LFC-1.names"

output_basename <- "../Outfiles/Avas_downDEGs_FDR_0.05_LFC-1_"#change accordingly 

#####
#
# read files
#
############

# list of DEGs
deg_names<-read.table(file=deg_namesPath)

##################
#
#
# Prep GO DB for DAG plotting
#
#
#
##################################################

#Key (col1) is children of (col2):
BP_Children <- as.data.frame(GOBPCHILDREN)

#BP_Children[BP_Children[,1] %in% GO_Term_Set[4],]

BP_Graph <- BP_Children[,c(2,1)]
colnames(BP_Graph)<-c("from","to")

#########################################
#
#
# GOTerm enrichment for Biological Process  
#
#
##################################################################

#all go terms
go_process<-readMappings(processPath)

module_transcripts_for_GO<-factor(as.integer(names(go_process) %in% deg_names$V1))
names(module_transcripts_for_GO)<-names(go_process)

table(module_transcripts_for_GO)

#for MF: ontology="MF"
#for BP: ontology="BP"
#for CC: ontology="CC"
process_go_data<-new("topGOdata", ontology="BP", allGenes=module_transcripts_for_GO, annot=annFUN.gene2GO, gene2GO=go_process, nodeSize=5)

process_resultsFisherClassic<-runTest(process_go_data, algorithm = "classic", statistic="fisher")

#for MF: topNodes=1182
#for BP: topNodes=6096
#for CC: topNodes=868
ProcessGoTable<-GenTable(process_go_data, classicFisher=process_resultsFisherClassic,
                         orderBy="classicFisher", ranksOf="classicFisher", topNodes=6096)


ProcessGoTable$classicFisher <- as.numeric(ProcessGoTable$classicFisher)

#check results
table(ProcessGoTable$classicFisher < 0.05)
table(ProcessGoTable$classicFisher < 0.01)
table(ProcessGoTable$classicFisher < 0.001)


#########################################
#
#
# Semantic similarity
#
#
##################################################################

selected_pvalue <- 0.05

BP_db <- godata('org.Avastus.eg.db', keytype = "GENENAME", ont = "BP")

GO_BP_SemanticSimilarity<- mgoSim(ProcessGoTable[ProcessGoTable$classicFisher < selected_pvalue,]$GO.ID, ProcessGoTable[ProcessGoTable$classicFisher < selected_pvalue,]$GO.ID, semData = BP_db, measure = "Wang", combine = NULL)

GO_BP_Tree <- hclust(dist(GO_BP_SemanticSimilarity, method = "euclidean"), method = "ward.D2")

#here you can change the minimum cluster size and the "granularity" increasing the deepSplit value to 2, 3, or 4 (4 favours smaller clusters).
#I find this setting is ok as it is
dynamic_cut <- cutreeDynamic(GO_BP_Tree, minClusterSize=5, method="hybrid", distM=as.matrix(dist(GO_BP_SemanticSimilarity, method = "euclidean")), deepSplit=1)
dynamic_colors<-labels2colors(dynamic_cut)
dynamic_cut_ordered <- dynamic_cut[order.dendrogram(as.dendrogram(GO_BP_Tree))] 

clusters_numbers <- unique(dynamic_cut_ordered) - (0 %in% dynamic_cut_ordered)
n_clusters <- length(clusters_numbers)

BP_dend <- branches_attr_by_clusters(as.dendrogram(GO_BP_Tree), dynamic_cut_ordered, values = standardColors(n=n_clusters)[unique(dynamic_cut_ordered)])
BP_dend <- sort(BP_dend)

#plot the heatmap, change name according to genes of interest!
pdf(file = paste0(output_basename,"goTerm_SematicSimilarity_Clusters.pdf"), paper = "a4", onefile = TRUE, width = 0, height = 0)

heatmap(GO_BP_SemanticSimilarity, Rowv = BP_dend, Colv = BP_dend)

dev.off()

####################
#
# prep a data.frame with all the results to save a csv
#
##################################################

BP_go_colors <- BP_dend %>% get_leaves_branches_col()

#plot colors
df<-as.data.frame(table(BP_go_colors))
pdf(paste0(output_basename,"cluster_color_key.pdf"), paper = "a4", onefile = TRUE, width = 0, height = 0,)#change the name according to GOIs
barplot(df$Freq~df$BP_go_colors,col=as.character(df$BP_go_colors),las=2,xlab = "")
dev.off()

BP_go_labels <- BP_dend %>% get_leaves_attr("label")

BP_go_groups <- data.frame(GO.ID = BP_go_labels, GO.colors = BP_go_colors, 
                           GO.Terms = select(GO.db, keys=BP_go_labels, columns = c("TERM")))

BP_go_groups <- merge(BP_go_groups, ProcessGoTable[ProcessGoTable$GO.ID %in% BP_go_groups$GO.ID,], by = "GO.ID")

#get ICs for the GO terms
GO_IC <- data.frame(GO.ID=names(BP_db@IC), GO.IC=BP_db@IC)
rownames(GO_IC) <- NULL

#
BP_go_groups <- merge(BP_go_groups, GO_IC[GO_IC$GO.ID %in% BP_go_groups$GO.ID,], by="GO.ID")

BP_go_groups$classicFisher <- as.numeric(BP_go_groups$classicFisher)


write.table(BP_go_groups, paste0(output_basename,"goTerm_SematicSimilarity_Clusters.csv"), dec = ",", sep = ";", row.names = FALSE)

####################
#
#
# Print chordPlot, DAG, and GO table for each GOTerm cluster
#
#
#
#################################################

pdf(file = paste0(output_basename,"_goTerm_SematicSimilarity_ClustersDetails.pdf"), paper = "a4", onefile = TRUE, width = 0, height = 0)

for (color in unique(BP_go_colors)) {
  
  go_title <- paste(color, "SemSim Cluster", color, sep = " ")
  
  GO_Term_Set <- BP_go_groups[BP_go_groups$GO.colors == color,]$GO.ID

  ################################################
  #
  #
  # DAG 
  #
  #
  #######################################################################
  
  BP_of_interest <- BP_Graph[BP_Graph$to %in% GO_Term_Set,]
  
  #set from to "" if the node is not in GO_Term_Set, i.e. is not enriched but part of the tree
  BP_of_interest$from <- ifelse(BP_of_interest$from %in% GO_Term_Set, BP_of_interest$from, "")
  
  #select only connected nodes
  BP_of_interest <- BP_of_interest[BP_of_interest$from != "",]
  
  #create graph with connected components
  connected_GOs <- unique(c(BP_of_interest$from,BP_of_interest$to))
  
  GO_Term_Set_Graph <- graph_from_data_frame(BP_of_interest, vertices = GO_IC[GO_IC$GO.ID %in% connected_GOs,])
  
  degree_zero_vertices <- names(degree(GO_Term_Set_Graph, mode = "in"))[degree(GO_Term_Set_Graph, mode = "in") == 0]
  
  write.table(degree_zero_vertices, paste0(output_basename,"_BP_",color,"_zeroDegreeVertices.txt"), dec = ".", sep = ";", row.names = FALSE, col.names = FALSE)
  
  
  #set vertex color
  grey_colors <- grey.colors(32*8,0,1)[seq(1,32*8,8)]
  
  vertex_IC <- vertex_attr(GO_Term_Set_Graph,"GO.IC")
  vertex_IC[which(vertex_IC == Inf)] <- max(vertex_IC[-which(vertex_IC == Inf)])
  vertex_IC <- (vertex_IC - min(vertex_IC)) / (max(vertex_IC)- min(vertex_IC))
  
  vertex.cols <- grey_colors[round(1+vertex_IC*31)]
  
  p_sizes <- -log(BP_go_groups[BP_go_groups$GO.ID %in% vertex_attr(GO_Term_Set_Graph,"name"),9])
  
  discretized_p_sizes <- round(8*(p_sizes - min(p_sizes))/(max(p_sizes)-min(p_sizes)),0)
  
  plot.igraph(simplify(GO_Term_Set_Graph), vertex.label.dist = 0.8, vertex.label.degree = -pi/2, vertex.color = vertex.cols, vertex.size = discretized_p_sizes, main = go_title, vertex.label.cex = 0.5, edge.arrow.size = 0.5)
  
  #############################
  #
  #
  # Print GO table
  #
  #connected_GOs
  ###############################################
  
  #selected_GOs <- BP_go_groups$GO.colors == color & BP_go_groups$GO.ID %in% connected_GOs
  
  #term_table <- BP_go_groups[selected_GOs,][,c(1,4,6:10)]
  #term_table[,2] <- sapply(apply(BP_go_groups[selected_GOs,][4], 1, strwrap, width = 80), paste, collapse="\n")
  #rownames(term_table) <- term_table[,1]
  #term_table[,1] <- NULL
  
  #t <- tableGrob(term_table, rows = rownames(term_table), cols = colnames(term_table), theme=ttheme_default(base_size = 4))
  
  #title_grob <- textGrob(go_title)
  
  #t <- gtable_add_rows(t, heights = grobHeight(title_grob) + unit(4,'mm'), pos = 0)
  #t <- gtable_add_grob(t, title_grob, 1, 1, 1, ncol(t), clip = "off")
  
  #grid.newpage()
  #grid.draw(t)
  
  #############
  #modified from https://jokergoo.github.io/circlize_book/book/advanced-usage-of-chorddiagram.html#customize-sector-labels
  #
  #rotated labels
  #
  ###################################
  
  
  
  GO_Term_Set_SemanticSimilarity <- as.matrix(GO_BP_SemanticSimilarity[rownames(GO_BP_SemanticSimilarity) %in% connected_GOs, colnames(GO_BP_SemanticSimilarity) %in% connected_GOs])  
  
  chordDiagram(100*GO_Term_Set_SemanticSimilarity, keep.diagonal = FALSE, symmetric = TRUE,
               link.visible = (GO_Term_Set_SemanticSimilarity >= 0.5), link.sort = TRUE, link.decreasing = FALSE,
               annotationTrack = c("grid","axis"), 
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(GO_Term_Set_SemanticSimilarity))))))
  
  # we go back to the first track and customize sector labels
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(-0.2, 0.5))
  }, bg.border = NA) # here set bg.border to NA is important
  
  title(go_title)
  
}

dev.off()

