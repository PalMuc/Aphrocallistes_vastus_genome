# synteny_2d_plot.R
# make dot plot of synteny between two genomes, based on unidirectional blast hits (i.e. not reciprocal)
# created by WRF 2019-04-01
# last modified 2022-11-07

#args = commandArgs(trailingOnly=TRUE)

# read data file from scaffold_synteny.py
#all2Dfile = args[1]
all2Dfile = "~/genomes/aphrocallistes_vastus_PORI/Aphrocallistes_vastus_genome/synteny/vs_oopsacas/Avas.1.29_vs_oopsacas_gb.scaffold2d_points.local.tab"

# read optional species names
# should match between query and subject in scaffold_synteny.py, as -f and -F

genome1_lab = "A. vastus (total Mb)"
genome2_lab = "O. minuta (total Mb)"


# read all data in a single file
all2Ddata = read.table(all2Dfile, sep="\t", stringsAsFactors=FALSE)
#head(all2Ddata)

# breakdown categories of either two genomes s1 and s2, or gene hits g
categories = all2Ddata[,1]

is_scaf1 = which(categories=="s1")
scafdata1 = all2Ddata[is_scaf1,]
#scafdata1

longestscaf1 = max(scafdata1[,6])
#longestscaf1
is_longscafs1 = which(as.numeric(scafdata1[,5]) > 0.0009)
#is_longscafs1

longscafs1 = c(0, scafdata1[,6][is_longscafs1] )
longscafs1_names = scafdata1[is_longscafs1,2]

is_scaf2 = which(categories=="s2")
scafdata2 = all2Ddata[is_scaf2,]
longestscaf2 = max(scafdata2[,6])
is_longscafs2 =  which(as.numeric(scafdata2[,5]) > 0.0009)
scafdata2_long = scafdata2[is_longscafs2,]
longscafs2 = c(0, scafdata2[,6][is_longscafs2] )
length(longscafs2)
longscafs2_names = scafdata2_long[,2]

is_points = which(categories=="g")
pointsdata = all2Ddata[is_points,]
is_both_longscaf = !is.na(match(pointsdata[,3], longscafs1_names)) & !is.na(match(pointsdata[,5], longscafs2_names))
pointsdata_long = pointsdata[is_both_longscaf,]
head(pointsdata_long)

scaffold_match_data = data.frame( sc1 = pointsdata_long[,3], sc2 = pointsdata_long[,5] )
match_freq_table = table(scaffold_match_data)

scaffold_best_match = apply(match_freq_table, 2, which.max)
scaffold_best_match
# moves alphabetic to length order
sc2_alpha_to_by_length = match( names(scaffold_best_match) , scafdata2[is_longscafs2,2] )

# moves length to alphabetic
sc2_by_length_to_alpha = match( scafdata2[is_longscafs2,2] , names(scaffold_best_match) )

sc2_sorted_index = sort(scaffold_best_match[sc2_by_length_to_alpha], index.return = TRUE)
scafdata2_reorder = scafdata2_long[sc2_sorted_index$ix,]
sc2_reorder_index = match( scafdata2_reorder[,2], names(sc2_sorted_index$x) )

# recount total length
longscafs1 = cumsum(c(0,as.numeric(scafdata1[is_longscafs1,4])))
longscafs2 = cumsum(c(0,as.numeric(scafdata2_reorder[,4])))
# reassign dots by scaffold
adj_points1_positions = pointsdata_long[,6] + longscafs1[match( pointsdata_long[,3], longscafs1_names ) ]
adj_points2_positions = pointsdata_long[,7] + longscafs2[match( pointsdata_long[,5], scafdata2_reorder$V2 ) ]


genome_x = adj_points2_positions
genome_y = adj_points1_positions
xmax = tail( pretty(longscafs2), n=1)
ymax = tail( pretty(longscafs1), n=1)
longscafs_x = longscafs2
longscafs_y = longscafs1
nscafs_x = length(longscafs2)
nscafs_y = length(longscafs1)
xlab = genome2_lab
ylab = genome1_lab


xmax_mb = round(xmax / 1000000)
ymax_mb = round(ymax / 1000000)

# larger bitscores make larger points
pointsize = log10(as.numeric(pointsdata[,8])) / 4

dotcolor = "#18519388" # color blue

# make PDF
#outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",all2Dfile,perl=TRUE)
outputfile = "~/genomes/aphrocallistes_vastus_PORI/Aphrocallistes_vastus_genome/synteny/vs_oopsacas/Avas.1.29_vs_oopsacas_gb.scaffold2d_points.local.small.pdf"
pdf(file=outputfile, width=5, height=7) # a4 size

par( mar=c(4.5,4.5,1,1) )

plot(genome_x, genome_y, pch=16, 
     xlim = c(0,60000000), ylim = c(0,80000000),
     col=dotcolor, cex=0.5, cex.lab=1.4, 
     main="", xlab=xlab, ylab=ylab, 
     axes=FALSE )

tickpoints = pretty(c(0,xmax_mb))
axis(1, at=tickpoints*1000000, labels=tickpoints, cex.axis=1.3)
barpos_x = rep(c( ymax*-0.01, ymax*-0.02),round(nscafs_x)/2)
#segments( longscafs_x[1:(nscafs_x-1)], barpos_x[1:(nscafs_x-1)], longscafs_x[2:nscafs_x], barpos_x[0:(nscafs_x-1)], lwd=3, col = "#88888888")
segments(longscafs_x, 0, longscafs_x, longscafs_y[nscafs_y], lwd=0.1, col="#88888888")

tickpoints = pretty(c(0,ymax_mb))
axis(2, at=tickpoints*1000000, labels=tickpoints, cex.axis=1.3, line=0.5)
barpos_y = rep(c( xmax*-0.01, xmax*-0.02),round(nscafs_y)/2)
#segments( barpos_y[1:(nscafs_y-1)], longscafs_y[1:(nscafs_y-1)], barpos_y[0:(nscafs_y-1)], longscafs_y[2:nscafs_y], lwd=3, col = "#88888888")
segments( 0, longscafs_y, longscafs_x[nscafs_x], longscafs_y, lwd=0.8, col="#00000088")

# make shortened names
short_names1 = gsub("Aphrocallistes_vastus_HiC-scaffold_", "", longscafs1_names)
short_names2 = gsub("JAKMXF010000", "", scafdata2_reorder[,2])

# display numbers beside axis scaffold segments
# this controls bars along the left side, so actually Y axis
max_x_n = 25 # max numbers to display
textpos_x = rep(c( ymax*-0.02, ymax*-0.01 ),round(max_x_n)/2)
textmidbar = as.numeric(scafdata1[1:max_x_n,6]) - as.numeric(scafdata1[1:max_x_n,4])/2
#text(textpos_x, textmidbar, short_names1[1:max_x_n], cex=0.5)

# this controls bars along the bottom side, so actually X axis
max_y_n = length(longscafs_x)-1
scaf2_over_1M = as.numeric(scafdata2_reorder[,4]) > 1000000
textpos_y = rep(c( xmax*-0.02, xmax*-0.01 ),round(max_y_n)/2)
textmidbar = ( as.numeric(longscafs_x[1:(max_y_n)]) + as.numeric(longscafs_x[2:(max_y_n+1)]) ) / 2
#text(textmidbar[scaf2_over_1M], textpos_y[scaf2_over_1M], short_names2[scaf2_over_1M], cex=0.5)

dev.off()

#
