# analysis of differential gene expression #
to re-run the DGE pipeline:
 
open the script `Avas_DeSeq2_Analysis.R`

make sure that the working directory is set to the directory `R/` (where the above mentioned script is stored)

source the script `Avas_DeSeq2_Analysis.R`

this will generate output files necessary for other R scripts.

after the script completes source the plot generating scripts: `Gen*Plots`

if more genelist are to be checked: modify the script `homolog_paralog_DEGs_makeLists.R` to include/remove the desired list in the graphs.

