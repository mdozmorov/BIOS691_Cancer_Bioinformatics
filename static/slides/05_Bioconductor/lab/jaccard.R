# Bioconducting the pairwise Jaccard workflow
# Define a function for the Jaccard statistic

jaccard <- function(x, y) {
gr_x <- import(x)
gr_y <- import(y)
intersects <- intersect(gr_x, gr_y, ignore.strand=TRUE)
unions <- union(gr_x, gr_y, ignore.strand=TRUE)
sum(width(intersects)) / sum(width(unions))
}

# Bioconducting the pairwise Jaccard workflow
# Compute the statistics in parallel

files <- Sys.glob("*.merge.bed")
jaccard_matrix <- outer(files, files,
function(a, b) mcmapply(jaccard, a, b))

# Make the plot
library(gplots)
library(RColorBrewer)
heatmap.2(jaccard_matrix, col = brewer.pal(9, "Blues"))
