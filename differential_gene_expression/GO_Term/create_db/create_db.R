# if installation is required, uncomment the following block
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("AnnotationForge")
#BiocManager::install("GO.db")
##

library("AnnotationForge")
library("GO.db")

#read
Genes <- read.csv("./avas_GOs4OwnDB_Genes.csv", sep = ";")
Chromo <- read.csv("./avas_GOs4OwnDB_Chromo.csv", sep = ";")
GOs <- read.csv("./avas_GOs4OwnDB_GOs.csv", sep = ";")

head(Genes)
head(Chromo)
head(GOs)

#check
table(GOs$GO %in% keys(GO.db))

#filter
Exclude_GIDs <- GOs[!(GOs$GO %in% keys(GO.db)),]$GID

f_Chromo <- Chromo[!(Chromo$GID %in% Exclude_GIDs),]
f_Genes <- Genes[!(Genes$GID %in% Exclude_GIDs),]
f_GOs <- GOs[!(GOs$GID %in% Exclude_GIDs),]

#check
table(f_GOs$GO %in% keys(GO.db))

## Then call the function
makeOrgPackage(gene_info=f_Genes, chromosome=f_Chromo, go=f_GOs,
               version="0.1",
               author = "Sergio Vargas <sergio.vargas@lmu.de>",
               maintainer = "Sergio Vargas <sergio.vargas@lmu.de>",
               outputDir = "./",
               tax_id = "83887",
               genus = "Aphrocallistes",
               species = "vastus",
               goTable="go")


