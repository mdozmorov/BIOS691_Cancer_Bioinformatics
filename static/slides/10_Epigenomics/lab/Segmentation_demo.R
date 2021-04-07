# Example taken from the nullranges package
# https://github.com/nullranges/nullranges
# BiocManager::install("nullranges/nullranges")
library(nullranges)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)

# Get list of genes
edb <- EnsDb.Hsapiens.v86
g <- genes(edb, filter = AnnotationFilterList(GeneIdFilter("ENSG", "startsWith")))
x <- keepSeqlevels(g, as.character(seq_len(12)), pruning.mode = "coarse")

# Tile the genome into 1M basepair windows
query <- tileGenome(seqlengths(x)[seqnames(x)@values], tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
# Calculate gene density per window
counts <- countOverlaps(query, x)
# Plot the results, with a little noise
eps <- rnorm(length(counts), 0, .2)
plot(sqrt(counts) + eps)

# HMM
library(RcppHMM)
hmm <- initPHMM(n = 3) # Three states
hmm <- learnEM(hmm,
  counts,
  iter = 400,
  delta = 1e-5,
  print = TRUE
)
hmm
# Find states with the Viterbi algorithm
v <- as.integer(factor(viterbi(hmm, counts), levels = hmm$StateNames))
plot(sqrt(counts) + eps, col = v)
