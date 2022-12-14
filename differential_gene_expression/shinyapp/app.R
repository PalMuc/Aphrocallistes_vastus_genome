# some processing of DGE for Aphrocallistes vastus
# shinyapp created by WRF 2022-11-15
# last modified 2022-11-26

library(shiny)
library(ggplot2)
library(dplyr)
library(ape)
library(ggtree)


# CHANGE PATHS HERE | | |
# CHANGE PATHS HERE v v v

# SET THIS: this must point to the functional annotation table
dge_data_file = "~/git/Aphrocallistes_vastus_genome/differential_gene_expression/Data/avas_v1_29_annotations_prot_fr_blastp_k10.annotation_table_clusters.tsv.gz"
# SET THIS: must point to files Avas_DGE_quant.csv and Avas_DGE_quant.info respectively
dge_raw_counts_file = "~/git/Aphrocallistes_vastus_genome/differential_gene_expression/Data/Avas_DGE_quant.csv.gz"
dge_sample_info_file = "~/git/Aphrocallistes_vastus_genome/differential_gene_expression/Data/Avas_DGE_quant.info"
# SET THIS: must point to the FOLDER containing all fasttree trees from the ortholog clusters
fasttree_folder = "~/git/Aphrocallistes_vastus_genome/ortholog_clusters/fasttree/"

# CHANGE PATHS HERE ^ ^ ^
# CHANGE PATHS HERE | | |

################################################################################
# read in all data files

# read final DGE table
dge_data = read.table(dge_data_file, header=TRUE, sep="\t", quote = '"')
#head(dge_data)

dge_raw_counts = read.table(dge_raw_counts_file, header = TRUE, sep = "\t")
dge_sample_info = read.csv(dge_sample_info_file , stringsAsFactors = FALSE)
dge_sample_xpos = match( dge_sample_info$individum, sort(unique(dge_sample_info$individum) ) )

#names(dge_data)
pchoices = names(dge_data)[24:29]
#pchoices
#[1] "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"         "padj"  

# make various things related to displaying trees
species_labels = c("AQUE", "AVAS", "BFLO", "CELE", "EMUE", "HCAL", "HEX06.EUSP", "HEX07.FOCC", 
                   "HHON", "HSAP", "MBRE", "NVEC", "RESC", "SCIL", "TADH" )
species_colors = c("#4292c6", #"AQUE"
                   "#00abff", #"AVAS"
                   "#fd8d3c", #"BFLO"
                   "#bf855c", #"CELE"
                   "#238443", #"EMUE"
                   "#b770df", #"HCAL"
                   "#22d6ef", #"HEX06-EUSP"
                   "#22d6ef", #"HEX07-FOCC"
                   "darkgrey", #"HHON"
                   "#804200", #"HSAP"
                   "#7a0177", #"MBRE"
                   "#a6014c", #NVEC
                   "#a6014c", #RESC
                   "#08306b", #"SCIL"
                   "gray" ) #"TADH"


# nothing needs to be set below
################################################################################

# begin app interface
ui <- fluidPage(
  titlePanel(basename(dge_data_file), windowTitle = "Avas v1.29 DGE explorer"),
  fluidRow(
    column(5,
           plotOutput(outputId = "pvaluePlot",
                      width="100%", height="650px",
                      click = "plot_click",
                      brush = brushOpts(id = "pvalue_brush")
           )
    ), # end column
    column(3,
           selectInput("x_select", label="X-axis control of plot",
                        choices = pchoices,
                        selected = "log2FoldChange"),
           selectInput("y_select", label="Y-axis control of plot",
                        choices = pchoices, 
                        selected = "baseMean"),
           sliderInput(inputId = "pvalueFilter",
                       label = "p significance (log)",
                       min = -20, max = -1, value = -4
           ),
           textInput("userGene", label = "Highlight gene model with full ID, e.g. Avas.s013.g133.i1", 
                     value = ""),
           textInput("userCluster", label = "Show members of cluster by name, e.g. paralog_cluster_04246", 
                     value = ""),
           textOutput("treeStats"),
           actionButton("treeUpdateButton", label = "Refresh cluster tree"),
           plotOutput(outputId = "rawScorePlot",
                      width="100%", height="200px" )
    ), # end column
    column(4,
           plotOutput(outputId = "selectedTree",
                      width="100%", height="700px"
           )
    ) # end column
  ),
  fluidRow(
    column(4,
           verbatimTextOutput("selectCountText"),
    ),
    tableOutput("selectedGenes")
  ) # end row
) # end fluidPage

server <- function(input, output) {
  #output$selectCountText <- renderPrint({ 10^(input$pvalueFilter) })
  #output$selectCountText <- renderPrint({ str(input$pvalue_brush) })
  
  output$pvaluePlot <- renderPlot({
    is_sig_point = dge_data$padj < 10^(input$pvalueFilter)
    #is_sig_point = dge_data$padj < 0.001

    # check if gene is in selected cluster
    is_user_cluster = ( dge_data$ortholog_cluster_id == input$userCluster )
    du = filter(dge_data, is_user_cluster)
    # take user defined genes, split into list, and check if row is in that list
    is_user_gene = ( dge_data$SeqName %in% as.list(scan(text=trimws( input$userGene ), what='', sep=','))  )
    dg = filter(dge_data, is_user_gene)
    
    point_color = ifelse(is_sig_point, "#045a8d88", "#fe992922")
    point_color[is.na(is_sig_point)] = "#00000022"
    
    gg = ggplot(data = dge_data, aes( x=get(input$x_select), y=get(input$y_select) )) +
      theme(legend.position="none",
            axis.text=element_text(size=16),
            axis.title=element_text(size=18),
            plot.title = element_text(size=23) ) +
      labs(x=input$x_select, y=input$y_select ) +
      scale_x_continuous(limits=c(-10,10)) +
      scale_y_log10() +
      #scale_color_manual(values=c("#fe992922","#045a8daa")) +
      #geom_point( aes(color=is_sig_point), size=point_size, shape=point_type) +
      geom_point( color=point_color, size=2 ) +
      geom_point(data = du, aes( x=get(input$x_select), y=get(input$y_select)),
                 color ="#dd0011", size=3, shape=5, stroke=4) + # show cluster members
      geom_point(data = dg, aes( x=get(input$x_select), y=get(input$y_select)),
                 color ="#0898ac", size=7, shape=16, alpha=0.6) + # show gene model
      annotate(geom="text", x=-Inf, y=Inf, size=5, hjust=0, vjust=1, 
               label=paste(table(is_sig_point)[["TRUE"]], "significant out of", length(is_sig_point),
                           "\nselected cluster contains", sum(is_user_cluster,na.rm = TRUE) ) )
      gg
  })
  
  tf <- eventReactive(input$treeUpdateButton, {
    # fix file names to match table
    full_cluster_name = gsub("paralog_cluster_", "paralogs/paralogs_holozoa_allprots_v10_", 
      gsub("homolog_cluster_" , "homologs/homologs_holozoa_allprots_v10_", input$userCluster) )
    tree_file = paste0(fasttree_folder, full_cluster_name, ".fasta.aln.tre" )
    read.tree(file = tree_file )
  })
  output$treeStats <- renderText({
    # get species counts for tree by parsing species names from seq IDs
    t = tf()
    tree_taxa = gsub("(^[\\w-]+)_(.*)","\\1",t$tip.label,perl=TRUE)
    paste( "Tree has", length( t$tip.label ), "sequences, for", length(unique(tree_taxa)), "species")
  })
  output$selectedTree <- renderPlot({
    t = tf()
    tree_taxa = gsub("(^[A-Za-z0-9-]+)_(.*)","\\1",t$tip.label,perl=TRUE)
    tree_cols = species_colors[match(tree_taxa, gsub("\\.","-",species_labels) )]
    ggtree( t ) +
      geom_tippoint(shape = 16, color=tree_cols, size=3, alpha=1 ) +
      geom_tiplab(color="black", show.legend = FALSE, size = 2.5,offset=0.04, angle=0) +
      geom_nodelab(size=2, nudge_x =0.04)
  })
  
  output$rawScorePlot <- renderPlot({
    # draw scatterplot showing the counts for each sample
    # body samples are left, tip samples are right
    gene_index = match( trimws(input$userGene), dge_raw_counts[,1])
    #gene_index = match( as.list(scan(text=trimws(input$userGene), what='', sep=','))[[1]] , dge_raw_counts[,1])
    gene_index = ifelse(!is.na(gene_index), gene_index, 1)
    dge_one_gene_df = data.frame( sample=dge_sample_xpos, counts=as.numeric(dge_raw_counts[gene_index,2:25]), condition=dge_sample_info$condition )
    ggplot(data = dge_one_gene_df, aes(x=sample, y=counts, color=condition)) +
      theme(legend.position="none",
            axis.title.x = element_blank(), axis.title.y = element_blank()) +
      scale_color_manual(values=c("#045a8d","#910f7c")) +
      scale_x_discrete(limits=c(1:6), labels=sort(unique(dge_sample_info$individum))) +
      scale_y_continuous(limits=c(0,max(dge_one_gene_df$counts))) +
      geom_point(size=5, alpha=0.7) +
      annotate(geom="text", x=2, y=Inf, label="body", vjust=1) +
      annotate(geom="text", x=5, y=Inf, label="tip", vjust=1)
  })
  
  output$selectedGenes <- renderTable({
    # combine user genes, brushed points, and user cluster into one table
    # this table can be redundant, meaning selected and brushed rows can repeat
    du = rbind( filter(dge_data, SeqName %in% as.list(scan(text=trimws( input$userGene ), what='', sep=',')) ),
               brushedPoints(dge_data, input$pvalue_brush, xvar = input$x_select, yvar = input$y_select ),
               filter(dge_data, ortholog_cluster_id == input$userCluster) )
  })

}

# Create Shiny app ----
shinyApp(ui = ui, server = server)

#