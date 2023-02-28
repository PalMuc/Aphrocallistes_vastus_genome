library(plotly)

balls=read.table(file="Aphrocallistes_reference_assembly_v1.ave.cov.gc.combined_no-mt.txt", header = TRUE)



# table headers:
#Name_of_contig
#scaffold_length
#PE_DNA_bowtie2
#PE_DNA_hisat2
#PE_RNA_hisat2
#ONT_RNA_gmap
#ONT_DNA_minimap2
#ONT_RNA_minimap2
#GC_percentage




## overview:


## DNA (short-reads) / RNA (short-reads)
# PE_DNA_bowtie2/ PE_RNA_hisat2 / GC_percentage:
plot_ly(balls, x = ~PE_DNA_bowtie2, y = ~PE_RNA_hisat2, z = ~GC_percentage, marker = list(color = ~PE_DNA_bowtie2, size = ~scaffold_length, sizeref= 1000, sizemode = "area", opacity=1, showscale = TRUE), text = ~paste('Scaffold name:', balls$Name_of_contig)) %>% add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PE_DNA_bowtie2',range = c(0,700), zerolinewidth = 2, ticklen = 1, gridwith = 1, zeroline = TRUE),
                      yaxis = list(title = 'PE_RNA_hisat2', range = c(0,2100), zerolinewidth = 2, ticklen = 1,gridwith = 1, zeroline = TRUE), 
                      zaxis = list(title = 'GC_percentage', range = c(0,70), zerolinewidth = 2, ticklen = 1, gridwith = 1, zeroline = TRUE), 
                      paper_bgcolor = "#444", plot_bgcolor = "#444"), annotations = list(x = 1.0, y = 1.0, text = 'PE_DNA_bowtie2avg. Coverage', showarrow = FALSE))  

## DNA (long-reads) / RNA (short-reads)
# ONT_DNA_minimap2/ PE_RNA_hisat2 / GC_percentage:
plot_ly(balls, x = ~ONT_DNA_minimap2, y = ~PE_RNA_hisat2, z = ~GC_percentage, marker = list(color = ~PE_DNA_bowtie2, size = ~scaffold_length, sizeref= 1000, sizemode = "area", opacity=1, showscale = TRUE), text = ~paste('Scaffold name:', balls$Name_of_contig)) %>% add_markers() %>%
  layout(scene = list(xaxis = list(title = 'ONT_DNA_minimap2', range = c(0,1000), zerolinewidth = 2, ticklen = 1, gridwith = 1, zeroline = TRUE),
                      yaxis = list(title = 'PE_RNA_hisat2', range = c(0,2100), zerolinewidth = 2, ticklen = 1,gridwith = 1, zeroline = TRUE), 
                      zaxis = list(title = 'GC_percentage', range = c(0,70), zerolinewidth = 2, ticklen = 1, gridwith = 1, zeroline = TRUE), 
                      paper_bgcolor = "#444", plot_bgcolor = "#444"), annotations = list(x = 1.0, y = 1.0, text = 'PE_DNA_bowtie2avg. Coverage', showarrow = FALSE))  

## DNA (short-reads) / RNA (long-reads)
# PE_DNA_bowtie2/ ONT_RNA_gmap / GC_percentage:
plot_ly(balls, x = ~PE_DNA_bowtie2, y = ~ONT_RNA_gmap, z = ~GC_percentage, marker = list(color = ~PE_DNA_bowtie2, size = ~scaffold_length, sizeref= 1000, sizemode = "area", opacity=1, showscale = TRUE), text = ~paste('Scaffold name:', balls$Name_of_contig)) %>% add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PE_DNA_bowtie2',range = c(0,700), zerolinewidth = 2, ticklen = 1, gridwith = 1, zeroline = TRUE),
                      yaxis = list(title = 'ONT_RNA_gmap', range = c(0,300), zerolinewidth = 2, ticklen = 1,gridwith = 1, zeroline = TRUE), 
                      zaxis = list(title = 'GC_percentage', range = c(0,70), zerolinewidth = 2, ticklen = 1, gridwith = 1, zeroline = TRUE), 
                      paper_bgcolor = "#444", plot_bgcolor = "#444"), annotations = list(x = 1.0, y = 1.0, text = 'PE_DNA_bowtie2avg. Coverage', showarrow = FALSE))  



# PE_DNA_bowtie2/ ONT_RNA_minimap2 / GC_percentage:
plot_ly(balls, x = ~PE_DNA_bowtie2, y = ~ONT_RNA_minimap2, z = ~GC_percentage, marker = list(color = ~PE_DNA_bowtie2, size = ~scaffold_length, sizeref= 1000, sizemode = "area", opacity=1, showscale = TRUE), text = ~paste('Scaffold name:', balls$Name_of_contig)) %>% add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PE_DNA_bowtie2',range = c(0,700), zerolinewidth = 2, ticklen = 1, gridwith = 1, zeroline = TRUE),
                      yaxis = list(title = 'ONT_RNA_minimap2', range = c(0,300), zerolinewidth = 2, ticklen = 1,gridwith = 1, zeroline = TRUE), 
                      zaxis = list(title = 'GC_percentage', range = c(0,70), zerolinewidth = 2, ticklen = 1, gridwith = 1, zeroline = TRUE), 
                      paper_bgcolor = "#444", plot_bgcolor = "#444"), annotations = list(x = 1.0, y = 1.0, text = 'PE_DNA_bowtie2avg. Coverage', showarrow = FALSE))  

## DNA (long-reads) / RNA (long-reads)
# ONT_DNA_minimap2/ ONT_RNA_gmap / GC_percentage:
plot_ly(balls, x = ~ONT_DNA_minimap2, y = ~ONT_RNA_gmap, z = ~GC_percentage, marker = list(color = ~PE_DNA_bowtie2, size = ~scaffold_length, sizeref= 1000, sizemode = "area", opacity=1, showscale = TRUE), text = ~paste('Scaffold name:', balls$Name_of_contig)) %>% add_markers() %>%
  layout(scene = list(xaxis = list(title = 'ONT_DNA_minimap2', range = c(0,1000), zerolinewidth = 2, ticklen = 1, gridwith = 1, zeroline = TRUE),
                      yaxis = list(title = 'ONT_RNA_gmap', range = c(0,300), zerolinewidth = 2, ticklen = 1,gridwith = 1, zeroline = TRUE), 
                      zaxis = list(title = 'GC_percentage', range = c(0,70), zerolinewidth = 2, ticklen = 1, gridwith = 1, zeroline = TRUE), 
                      paper_bgcolor = "#444", plot_bgcolor = "#444"), annotations = list(x = 1.0, y = 1.0, text = 'PE_DNA_bowtie2avg. Coverage', showarrow = FALSE))  

## DNA (long-reads) / RNA (long-reads)
# ONT_DNA_minimap2/ ONT_RNA_minimap2 / GC_percentage:
plot_ly(balls, x = ~ONT_DNA_minimap2, y = ~ONT_RNA_minimap2, z = ~GC_percentage, marker = list(color = ~PE_DNA_bowtie2, size = ~scaffold_length, sizeref= 1000, sizemode = "area", opacity=1, showscale = TRUE), text = ~paste('Scaffold name:', balls$Name_of_contig)) %>% add_markers() %>%
  layout(scene = list(xaxis = list(title = 'ONT_DNA_minimap2', range = c(0,1000), zerolinewidth = 2, ticklen = 1, gridwith = 1, zeroline = TRUE),
                      yaxis = list(title = 'ONT_RNA_minimap2', range = c(0,300), zerolinewidth = 2, ticklen = 1,gridwith = 1, zeroline = TRUE), 
                      zaxis = list(title = 'GC_percentage', range = c(0,70), zerolinewidth = 2, ticklen = 1, gridwith = 1, zeroline = TRUE), 
                      paper_bgcolor = "#444", plot_bgcolor = "#444"), annotations = list(x = 1.0, y = 1.0, text = 'PE_DNA_bowtie2avg. Coverage', showarrow = FALSE))  





## RNA (short-reads) / RNA (long-reads)
# PE_RNA_hisat2/ ONT_RNA_gmap / GC_percentage:
plot_ly(balls, x = ~PE_RNA_hisat2, y = ~ONT_RNA_gmap, z = ~GC_percentage, marker = list(color = ~PE_DNA_bowtie2, size = ~scaffold_length, sizeref= 1000, sizemode = "area", opacity=1, showscale = TRUE), text = ~paste('Scaffold name:', balls$Name_of_contig)) %>% add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PE_RNA_hisat2',range = c(0,500), zerolinewidth = 2, ticklen = 1, gridwith = 1, zeroline = TRUE),
                      yaxis = list(title = 'ONT_RNA_gmap', range = c(0,300), zerolinewidth = 2, ticklen = 1,gridwith = 1, zeroline = TRUE), 
                      zaxis = list(title = 'GC_percentage', range = c(0,70), zerolinewidth = 2, ticklen = 1, gridwith = 1, zeroline = TRUE), 
                      paper_bgcolor = "#444", plot_bgcolor = "#444"), annotations = list(x = 1.0, y = 1.0, text = 'PE_DNA_bowtie2avg. Coverage', showarrow = FALSE))  

# PE_RNA_hisat2/ ONT_RNA_minimap2 / GC_percentage:
plot_ly(balls, x = ~PE_RNA_hisat2, y = ~ONT_RNA_minimap2, z = ~GC_percentage, marker = list(color = ~PE_DNA_bowtie2, size = ~scaffold_length, sizeref= 1000, sizemode = "area", opacity=1, showscale = TRUE), text = ~paste('Scaffold name:', balls$Name_of_contig)) %>% add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PE_RNA_hisat2',range = c(0,500), zerolinewidth = 2, ticklen = 1, gridwith = 1, zeroline = TRUE),
                      yaxis = list(title = 'ONT_RNA_minimap2', range = c(0,300), zerolinewidth = 2, ticklen = 1,gridwith = 1, zeroline = TRUE), 
                      zaxis = list(title = 'GC_percentage', range = c(0,70), zerolinewidth = 2, ticklen = 1, gridwith = 1, zeroline = TRUE), 
                      paper_bgcolor = "#444", plot_bgcolor = "#444"), annotations = list(x = 1.0, y = 1.0, text = 'PE_DNA_bowtie2avg. Coverage', showarrow = FALSE))  



## RNA (long-reads) / RNA (long-reads)
# ONT_RNA_minimap2/ ONT_RNA_gmap / GC_percentage:
plot_ly(balls, x = ~ONT_RNA_minimap2, y = ~ONT_RNA_gmap, z = ~GC_percentage, marker = list(color = ~PE_DNA_bowtie2, size = ~scaffold_length, sizeref= 1000, sizemode = "area", opacity=1, showscale = TRUE), text = ~paste('Scaffold name:', balls$Name_of_contig)) %>% add_markers() %>%
  layout(scene = list(xaxis = list(title = 'ONT_RNA_minimap2', range = c(0,300), zerolinewidth = 2, ticklen = 1, gridwith = 1, zeroline = TRUE),
                      yaxis = list(title = 'ONT_RNA_gmap', range = c(0,300), zerolinewidth = 2, ticklen = 1,gridwith = 1, zeroline = TRUE), 
                      zaxis = list(title = 'GC_percentage', range = c(0,70), zerolinewidth = 2, ticklen = 1, gridwith = 1, zeroline = TRUE), 
                      paper_bgcolor = "#444", plot_bgcolor = "#444"), annotations = list(x = 1.0, y = 1.0, text = 'PE_DNA_bowtie2avg. Coverage', showarrow = FALSE))  


