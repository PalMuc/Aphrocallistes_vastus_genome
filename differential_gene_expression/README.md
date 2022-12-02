# differential gene expression #
Reads were mapped to transcripts with [salmon](https://combine-lab.github.io/salmon/), by [Patro et al 2017](https://doi.org/10.1038/nmeth.4197)

### analysis pipeline ###
to re-run the DGE pipeline:
 
open the script `Avas_DeSeq2_Analysis.R`

make sure that the working directory is set to the directory `R/` (where the above mentioned script is stored)

source the script `Avas_DeSeq2_Analysis.R`

this will generate output files necessary for other R scripts.

after the script completes source the plot generating scripts: `Gen*Plots`

if more genelist are to be checked: modify the script `homolog_paralog_DEGs_makeLists.R` to include/remove the desired list in the graphs.

### interactive DGE viewer ###
This is an app using [Rshiny](https://shiny.rstudio.com/), which can be run in Rstudio. The goal is to interactively view and select significantly up- or down-regulated points or clusters, to identify any patterns among the thousands of points. For instance, it would be useful to highlight all members of an ortholog cluster to see if some members of a cluster are upregulated, and others are downregulated. This can, in theory, be used by anyone, or modified for your own clusters or genomes. Please contact [user WRF](https://github.com/wrf) with questions about implementation.

![dge_viewer_shinyapp_screenshot_v1.png](https://github.com/PalMuc/Aphrocallistes_vastus_genome/blob/main/differential_gene_expression/shinyapp/dge_viewer_shinyapp_screenshot_v1.png)
