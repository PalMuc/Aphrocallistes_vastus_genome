# some processing of DGE for Spongilla lactustris single cell
# data from Musser 2021 Science
# shinyapp created by WRF 2022-12-23
# last modified 2022-12-23

library(shiny)
library(ggplot2)
library(dplyr)

################################################################################

spongilla_data_file = "~/project/musser2021_sponge-single-cell/supplementary_files/Data_S1_sheet3_diff_exp_cell_types.tab"
spongilla_data = read.table(spongilla_data_file, header=TRUE, sep="\t", quote='"', 
                            stringsAsFactors = FALSE, 
                            colClasses = c("character","character","character","character","numeric","numeric","numeric"))
pchoices = names(spongilla_data)[4:7]
cellcounts = table(spongilla_data$Cell.Type)
celltypes = names(cellcounts)


################################################################################

ui <- fluidPage(
  titlePanel(basename(spongilla_data_file), windowTitle = "Spongilla single-cell RNAseq explorer"),
  fluidRow(
    column(3,
           selectInput("celltype", label="Select cell type for right plot",
                       choices = celltypes,
                       selected = "Scl"),
    ), # end column
    column(3,
           textInput("userGene", label = "Highlight gene model with full ID, e.g. c13958-g1", 
                     value = ""),
    ), # end column
    column(3,
           sliderInput(inputId = "xZoom",
                       label = "X-axis range",
                       min = 0, max = 1, value = c(0,5) )
    ), # end column
    column(3,
           sliderInput(inputId = "yZoom",
                       label = "Y-axis range",
                       min = 0, max = 1, value = c(0,1) )
    ) # end column
  ),
  fluidRow(
    column(6,
           plotOutput(outputId = "allCellPlot",
                      width="100%", height="600px",
                      brush = brushOpts(id = "allcell_brush")
           )
    ), # end column
    column(6,
           plotOutput(outputId = "chosenCellPlot",
                      width="100%", height="600px",
                      brush = brushOpts(id = "chosen_brush")
           )
    ) # end column
  ),
  fluidRow(
    column(4,
           verbatimTextOutput("selectText"),
    ),
    tableOutput("selectedGenes")
  ) # end row
) # end fluidPage


server <- function(input, output) {
  output$selectCountText <- renderPrint({ cellcounts[match(input$celltype,celltypes)] })
  
  output$allCellPlot <- renderPlot({
    # check if gene is in selected cell type
    is_user_celltype = ( spongilla_data$Cell.Type == input$celltype )
    du = filter(spongilla_data, is_user_celltype)
    
    gga = ggplot(spongilla_data, aes(x=Percent.cells.expressing..outside.cluster., y=Percent.cells.expressing..in.cluster.)) +
      theme(legend.position="none",
            axis.text=element_text(size=16),
            axis.title=element_text(size=18),
            plot.title = element_text(size=23) ) +
      scale_x_continuous(limits=input$xZoom, expand = c(0.02,0) ) +
      scale_y_continuous(limits=input$yZoom, expand = c(0.02,0) ) +
      geom_point(aes(size=Mean.Log.Fold.Diff. , color=Cell.Type), alpha=0.3) +
      geom_point(data = du, aes(x=Percent.cells.expressing..outside.cluster., y=Percent.cells.expressing..in.cluster.),
                 color ="#000000", size=2) # show cell type
    gga
  })
  
  output$chosenCellPlot <- renderPlot({
    # check if gene is in selected cell type
    is_user_celltype = ( spongilla_data$Cell.Type == input$celltype )
    du = filter(spongilla_data, is_user_celltype)
    
    ggc = ggplot(du, aes(x=Percent.cells.expressing..outside.cluster., y=Percent.cells.expressing..in.cluster.)) +
      theme(legend.position="none",
            axis.text=element_text(size=16),
            axis.title=element_text(size=18),
            plot.title = element_text(size=23) ) +
      scale_x_continuous(limits=input$xZoom, expand = c(0.02,0) ) +
      scale_y_continuous(limits=input$yZoom, expand = c(0.02,0) ) +
      scale_color_continuous(type="gradient", high="#dddddd", low="#67001f") +
      geom_point(aes(size=Mean.Log.Fold.Diff.), color="#67001f", alpha=0.5)
    ggc
  })
  
  output$selectedGenes <- renderTable({
    # combine brushed points from both plots into one table
    is_user_celltype = ( spongilla_data$Cell.Type == input$celltype )
    du = filter(spongilla_data, is_user_celltype)
    rbind( brushedPoints(du, input$chosen_brush, 
                              xvar = "Percent.cells.expressing..outside.cluster." , 
                              yvar = "Percent.cells.expressing..in.cluster." ),
                brushedPoints(spongilla_data, input$allcell_brush, 
                              xvar = "Percent.cells.expressing..outside.cluster." , 
                              yvar = "Percent.cells.expressing..in.cluster." )
               )
  })#, digits=-2)
  
}

shinyApp(ui = ui, server = server)

#