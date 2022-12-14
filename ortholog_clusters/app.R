# multi-way interactive plot of orthologs from all-v-all MCL clustering
#
# last modified 2022-10-25

library(shiny)
library(ggplot2)
library(dplyr)
library(ape)
library(ggtree)

#install.packages("BiocManager")
#install.packages("ggfun")
#remotes::install_github("YuLab-SMU/ggtree")

cluster_table_headers = c("clusterName", "num_sequences", "num_taxa",
                          "min_per_taxon", "med_per_taxon", "max_per_taxon",
                          "min_length", "mean_length", "max_length", "taxaList")

################################################################################

### change this to output from makehomologs.py ###
input_cluster_file = "~/git/Aphrocallistes_vastus_genome/ortholog_clusters/fasta_clusters.H.holozoa_allprots_v10.tab.gz"
### change this to the fasttree tre directory ###
fasttree_folder = "~/git/Aphrocallistes_vastus_genome/ortholog_clusters/fasttree/"
pfam_folder = "~/git/Aphrocallistes_vastus_genome/ortholog_clusters/pfam/"

# DEBUG
# domain_file = "~/git/pfam_v8/paralogs_E_holozoa_allprots_v8_00097.pfam.tab"
# domain_data = read.table(domain_file, sep="\t")
# 
# domain_names = as.character(domain_data$V1)
# protein_names = as.character(unique(domain_data$V4))
# domain_evalue = as.numeric(domain_data$V13)
# 
# protein_index = match(domain_data$V4,protein_names)
# domain_index = protein_index + match(domain_data$V1, unique(domain_names))/20
# 
# domain_starts = as.integer(domain_data$V18)
# domain_ends = as.integer(domain_data$V19)
# domain_middle = (domain_starts + domain_ends ) / 2
# domain_width = domain_ends - domain_starts
# block_height = rep(0.4, length(protein_index))
# 
# dom_df = filter( 
#     data.frame(protein_index, domain_middle, domain_index, domain_width, block_height, domain_names, domain_evalue),
#     protein_index <= domplot_max_sp &
#         domain_evalue < 0.01
# )
# 
# ggplot(dom_df, aes(x=domain_middle, y=domain_index , 
#                    width=domain_width, height=block_height ,
#                    fill=domain_names) ) +
#     theme(legend.position = "top") +
#     labs(x=NULL, y=NULL) +
#     geom_tile(alpha=0.9) +
#     annotate("text", x=rep(0, max(dom_df$protein_index)), y=seq(0.9,max(dom_df$protein_index),1), 
#              label=protein_names[1:max(dom_df$protein_index)], hjust="left" )

################################################################################


#cluster_data = read.table(input_cluster_file, header=FALSE, sep="\t", stringsAsFactors = FALSE,
#                          col.names = cluster_table_headers )
cluster_data = read.table(input_cluster_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
cluster_table_headers = colnames(cluster_data)

clust_max_col = ifelse(ncol(cluster_data)==10,9,ncol(cluster_data))
clust_max_col
species_labels = cluster_table_headers[10:ncol(cluster_data)]
species_labels
species_count = length(species_labels)
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


# for 10 sp "AQUE" "AVAS" "BFLO" "EMUE" "HHON" 
#           "HSAP" "MBRE" "SCIL" "SROS" "TADH"
#species_colors = c("#4292c6", "#00abff", "#fd8d3c", "#238443", "darkgrey",
#                   "#800026", "#7a0177", "#08306b", "#3f007d", "gray" )

#domplot_max_sp = 20 # cap at displaying 20 proteins with domains, to not clog view

d = mutate(cluster_data[,1:clust_max_col],
           average_seqs_per_sp = num_sequences/num_taxa,
           median_vs_min = (med_per_taxon - min_per_taxon),
           max_vs_median = (max_per_taxon - med_per_taxon),
           is_avg_one = (num_sequences==num_taxa)
           )

ui <- fluidPage(
    titlePanel("Ortholog clusters all-vs-all"),
    h4(input_cluster_file),
    fluidRow(
        # subcolumns for input selectors
        column(3,
               selectInput("x_select_left", label="X-axis control of upper plot",
                           choices = names(d)[2:(clust_max_col+3)],
                           selected = "num_sequences"),
               selectInput("y_select_left", label="Y-axis control of upper plot",
                           choices = names(d)[2:(clust_max_col+3)], 
                           selected = "num_taxa")
        ),
        column(3,
               selectInput("x_select_right", label="X-axis control of lower plot",
                           choices = names(d)[2:(clust_max_col+3)],
                           selected = "num_sequences"),
               selectInput("y_select_right", label="Y-axis control of lower plot",
                           choices = names(d)[2:(clust_max_col+3)], 
                           selected = "average_seqs_per_sp")
        ),
        column(3,
               textInput("userCluster", label = "Cluster name, excluding .tree or .fasta", 
                         value = "paralog_cluster_00107"),
               textOutput("treeStats"),
               br(),
               actionButton("treeUpdateButton", label = "Refresh tree and domains")
        ),
        column(3,
               # numericInput("startingProtIndex", label = "index of first protein", 
               #              value = 1, min=1),
               numericInput("domainMaxSeqs", label = "max proteins to show in domain plot", 
                            value = 20, min=1)
        )
    ),
    fluidRow(
        column(4,
            # plots
           plotOutput(outputId = "seqTaxaCounts",
                      width="100%",
                      brush = "seqTaxaBrush"
           ),
           plotOutput(outputId = "maxperTaxon",
                      width="100%",
                      brush = "maxseqsBrush"
           )
        ), # end column
        column(4,
               plotOutput(outputId = "selectedTree",
                          width="100%", height="800px"
               )
        ), # end column
        column(4,
               plotOutput(outputId = "domainPlot",
                          width="100%", height="750px"
               )
        ) # end column
    ), # end fluidRow
    fluidRow(
        column(4,
               textOutput("selectCountText"),
               ),
        tableOutput("selectedClusters")
    )
) # end fluidPage

# Define server logic ----
server <- function(input, output) {
    
    # setup for combining brushed points from all 4 plots
    vals = reactiveValues(
        cls = d[NULL,]
    )
    observeEvent( c(input$seqTaxaBrush, input$avgProtBrush, input$maxseqsBrush, input$MMMBrush), {
        vals$cls = rbind(
            brushedPoints(d, input$seqTaxaBrush, xvar = input$x_select_left, yvar = input$y_select_left ),
            #brushedPoints(d, input$avgProtBrush, xvar = "num_sequences", yvar = "average_seqs_per_sp"),
            #brushedPoints(d, input$maxseqsBrush, xvar = "num_sequences", yvar = "mean_length"),
            brushedPoints(d, input$maxseqsBrush, input$x_select_right, yvar = input$y_select_right )
        )
    })

    # display the 4 plots, and add black dots of the selected points
    output$seqTaxaCounts <- renderPlot({
        ggplot(d, aes(x=get(input$x_select_left), y=get(input$y_select_left), color=is_avg_one)) +
            labs(x=input$x_select_left, y=input$y_select_left ) +
            theme(legend.position = "none") +
            geom_point(size=3) +
            geom_point(data=vals$cls, colour="#000000")
    })
    # output$avgProtperSp <- renderPlot({
    #      ggplot(d, aes(x=num_sequences, y=average_seqs_per_sp, color=is_avg_one)) +
    #         theme(legend.position = "none") +
    #         geom_point(size=3) +
    #         geom_point(data=vals$cls, colour="#000000")
    # })
    output$maxperTaxon <- renderPlot({
        ggplot(d, aes(x=get(input$x_select_right), y=get(input$y_select_right), color=is_avg_one)) +
            labs(x=input$x_select_right, y=input$y_select_right ) +
            theme(legend.position = "none") +
            geom_point(size=3) +
            geom_point(data=vals$cls, colour="#000000")
    })

    tf <- eventReactive(input$treeUpdateButton, {
      full_cluster_name = gsub("paralog_cluster_", "paralogs/paralogs_holozoa_allprots_v10_", 
                               gsub("homolog_cluster_" , "homologs/homologs_holozoa_allprots_v10_", input$userCluster) )
      tree_file = paste0(fasttree_folder, full_cluster_name, ".fasta.aln.tre" )
      read.tree(file = tree_file )
    })
    output$treeStats <- renderText({
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
    
    df <- eventReactive(input$treeUpdateButton, {
      full_cluster_name = gsub("paralog_cluster_", "paralogs/paralogs_holozoa_allprots_v10_", 
                               gsub("homolog_cluster_" , "homologs/homologs_holozoa_allprots_v10_", input$userCluster) )
      domain_file = paste0(pfam_folder, full_cluster_name, ".pfam.tab" )
    })
    output$domainPlot <- renderPlot({
        domain_file = df()
        domain_data = read.table(domain_file, sep="\t")
        
        domain_names = as.character(domain_data$V1)
        protein_names = as.character(unique(domain_data$V4))
        domain_evalue = as.numeric(domain_data$V13)
        num_prots = length(protein_names)
        
        protein_index = match(domain_data$V4,protein_names)
        domain_index = protein_index + 0.4*match(domain_data$V1, unique(domain_names))/length(unique(domain_names))
        
        domain_starts = as.integer(domain_data$V18)
        domain_ends = as.integer(domain_data$V19)
        domain_middle = (domain_starts + domain_ends ) / 2
        domain_width = domain_ends - domain_starts
        block_height = rep(0.4, length(protein_index))
        
        dom_df_unf = data.frame(protein_index, domain_middle, domain_index, domain_width, block_height, domain_names, domain_evalue)
        
        max_seqs_for_display = ifelse(input$domainMaxSeqs > num_prots, num_prots, input$domainMaxSeqs)
        min_prot_index = 1
        max_prot_index = (1 + max_seqs_for_display - 1)
        dom_df = filter( dom_df_unf, 
            .data$protein_index <= max_prot_index &
            .data$protein_index >= min_prot_index &
                domain_evalue < 0.01
        )
        
        ggplot(dom_df, aes(x=domain_middle, y=domain_index , 
                           width=domain_width, height=block_height ,
                           fill=domain_names) ) +
            theme(legend.position = "top" ) +
            labs(x=NULL, y=NULL,
                 fill="Domain name:") +
            geom_tile(alpha=0.9) +
            annotate("text", x=rep(0, max_seqs_for_display), y=seq(0.9,max_seqs_for_display,1), 
                     label=protein_names[min_prot_index:max_prot_index], hjust="left" )

    })
    
    # print row count of selected points
    output$selectCountText <- renderText({
        paste(nrow(vals$cls), "points selected")
    })
    
    # make table of the selected points across all 4 plots
    output$selectedClusters <- renderTable({
        vals$cls
    })
}

shinyApp(ui = ui, server = server)

#