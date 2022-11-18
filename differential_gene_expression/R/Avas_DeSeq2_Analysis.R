#http://www-huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html#preparing-count-matrices consulted 07.03.2016

#####################################
# requires R-4.1.2
#
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#
#BiocManager::install(c("DESeq2","geneplotter", "BiocParallel", "gplots" ))
#install.packages(c("RcolorBrewer","ggplot2"))
#
#
#################################################

library("DESeq2")
library("geneplotter")
library("gplots")
library("RColorBrewer")
library("BiocParallel")
library("ggplot2")
register(MulticoreParam(8))

######################
#
# Private functions
# 
##############################################

#######
#
# Private variables
#
##############################


#ref level for DESeq comparisons
referenceLevel<-"Body"

#what to compare? this is dataset specific and based on your .info file
compareThis<-"condition"

#should plots be produced?
withPlots<-FALSE

#levels of the condition to be tested
conditionLevels<-c("Body", "Tip")

#padj cutoff to use
padjLevel<-0.05

#LFC cutoff to use
lfcLevel<-1

#write output files
write_files<-TRUE

#######
#
# Path to files to be read in
#
##############################
input_counts<-"../Data/Avas_DGE_quant.csv"
input_sample_information<-"../Data/Avas_DGE_quant.info"

########
#
# Paths to store output files and plots. NO / at the end of the path!!!!!!
#
####################

outPath<-"../Outfiles"

##############
#
# Main: do not change anything from here unless you know what you are doing
#
########################

Avas_DE<-read.csv(input_counts, head=T, sep="\t")
Avas_DE_INFO<-read.csv(input_sample_information, head=T)

rownames(Avas_DE)<-Avas_DE$Name
Avas_DE$Name<-NULL

# remove the non-transdecoded transcripts
head(Avas_DE)

#truncate pseudocounts
Avas_DE<-trunc(Avas_DE, 0)

#genes with zero counts over all samples
table(rowSums(Avas_DE) > 0)
# Remove genes that have no counts over all samples
Avas_DE <- Avas_DE[(rowSums(Avas_DE) > 0),]

#read the data into an DESeq2 object
Avas_DeSeq<-DESeqDataSetFromMatrix(countData=Avas_DE, colData=Avas_DE_INFO, design=~condition)

#Collapse tech reps
Avas_DeSeq <- collapseReplicates(Avas_DeSeq, Avas_DeSeq$individum)

#relevel the factors so that control is the reference
Avas_DeSeq$condition<-relevel(Avas_DeSeq$condition,ref=referenceLevel)

#estimate size factors
#Thus, if all size factors are roughly equal to one, the libraries have been sequenced equally deeply.
Avas_DeSeq<-estimateSizeFactors(Avas_DeSeq)

#Check size factors if wanted values should be near zero.
Avas_DeSeq$sizeFactor-1


#get rows with all non-zero counts
non_zero_rows<-apply(counts(Avas_DeSeq), 1, function(x){all(x>0)})

#get all rows
all_rows<-apply(counts(Avas_DeSeq), 1, function(x){all(x>=0)})

#number of non-zero rows, i.e. transcript was detected in all replicates:
sum(non_zero_rows)

#number of rows, i.e. total number of rows:
sum(all_rows)


if(withPlots){

  #cummulative distribution of normalized counts for non-zero rows
  multiecdf(counts(Avas_DeSeq, normalized=T)[non_zero_rows, ], xlab="mean counts", xlim=c(0,1000))
  
  #density of normalized counts
  multidensity(counts(Avas_DeSeq, normalized=T)[non_zero_rows, ], xlab="mean counts", xlim=c(0,1000))
  
  #compare samples pair-wise if wanted. Uncomment and modify as required
  #MDPlot(counts(Avas_DeSeq, normalized=T)[non_zero_rows, ],c(1,2), main=paste( colnames(Avas_DeSeq)[1], "vs.", colnames(Avas_DeSeq)[2]))
}

#transform the data to rlog
Avas_DeSeq_rlog<-rlogTransformation(Avas_DeSeq, blind = T)

#calculate transformed distances
Avas_distances<-dist(t(assay(Avas_DeSeq_rlog)))

if(withPlots){
  #produce a heat plot using the transformed distances
  heatmap.2(as.matrix(Avas_distances), trace="none", col=rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)))

  #PCA of samples
  plotPCA(Avas_DeSeq_rlog, intgroup=(compareThis))
}


###
#
#Differential expression analysis
#
#####

Avas_DeSeq$condition<-factor(Avas_DeSeq$condition, levels = conditionLevels)

Avas_DeSeq<-estimateSizeFactors(Avas_DeSeq)
Avas_DeSeq<-estimateDispersions(Avas_DeSeq)

#plotDispEsts(Avas_DeSeq)

#wald test
Avas_DeSeq<-nbinomWaldTest(Avas_DeSeq)
Avas_DeSeq_Results<-results(Avas_DeSeq, pAdjustMethod = "BH")

#global transcriptome
table(Avas_DeSeq_Results$padj < padjLevel)
table(Avas_DeSeq_Results$padj < padjLevel  & abs(Avas_DeSeq_Results$log2FoldChange) >= lfcLevel)

#plot of p values to see if there are any notworthy trends
#hist(Avas_DeSeq_Results$pvalue, main = "Body vs. Tip", xlab="p-values")

#plotMA(Avas_DeSeq_Results, alpha=0.01)

#get significant genes
Avas_DeSeq_Results_DEGs<-subset(Avas_DeSeq_Results, padj < padjLevel)
Avas_DeSeq_Results_DEGs_names<-rownames(Avas_DeSeq_Results_DEGs)

#upregulated
Avas_DeSeq_Results_upDEGs<-subset(Avas_DeSeq_Results, padj < padjLevel & log2FoldChange >= lfcLevel)
Avas_DeSeq_Results_upDEGs_names<-rownames(Avas_DeSeq_Results_upDEGs)

#downregulated
Avas_DeSeq_Results_downDEGs<-subset(Avas_DeSeq_Results, padj < padjLevel & log2FoldChange <= -lfcLevel)
Avas_DeSeq_Results_downDEGs_names<-rownames(Avas_DeSeq_Results_downDEGs)


if(write_files){

  write.csv(Avas_DeSeq_Results, paste(outPath,"/Avas_AllTranscripts_DESeq2Results.csv", sep=""))
  write.csv(Avas_DeSeq_Results_DEGs, paste(outPath,"/Avas_DEGs_FDR_0.05.csv", sep=""))
  write.csv(Avas_DeSeq_Results_DEGs_names, paste(outPath,"/Avas_DEGs_FDR_0.05.names", sep=""),row.names = F,col.names = F)
  
  #up
  write.csv(Avas_DeSeq_Results_upDEGs, paste(outPath,"/Avas_upDEGs_FDR_0.05_LFC1.csv", sep=""))
  write.csv(Avas_DeSeq_Results_upDEGs_names, paste(outPath,"/Avas_upDEGs_FDR_0.05_LFC1.names", sep=""),row.names = F,col.names = F)
  
  #down
  write.csv(Avas_DeSeq_Results_downDEGs, paste(outPath,"/Avas_downDEGs_FDR_0.05_LFC-1", sep=""))
  write.csv(Avas_DeSeq_Results_downDEGs_names, paste(outPath,"/Avas_downDEGs_FDR_0.05_LFC-1.names", sep=""),row.names = F,col.names = F)
    
  #counts
  write.csv(counts(Avas_DeSeq), paste(outPath,"/Avas_DESeq2_counts.csv", sep=""))
}


