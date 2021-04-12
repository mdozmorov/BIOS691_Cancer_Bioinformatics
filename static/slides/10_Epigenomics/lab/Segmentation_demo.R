# Example taken from the nullranges package
# https://github.com/nullranges/nullranges
# BiocManager::install("nullranges/nullranges")
library(RcppHMM) # For Hidden Markov Models
library(DNAcopy) # for Circular Binary Segmentation
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)

# Get list of genes
edb <- EnsDb.Hsapiens.v86
g <- genes(edb, filter = AnnotationFilterList(GeneIdFilter("ENSG", "startsWith")))
x <- keepSeqlevels(g, c("2"), pruning.mode = "coarse") # Chr2 only as example

# Tile the genome into windows
query <- tileGenome(seqlengths(x)[seqnames(x)@values], tilewidth = 5e5, cut.last.tile.in.chrom = TRUE)
# Calculate gene density per window
counts <- countOverlaps(query, x)
# Plot the results
plot(sqrt(counts))

# HMM
hmm <- initPHMM(n = 3) # Three states
hmm <- learnEM(hmm,
               counts,
               iter = 400,
               delta = 1e-5,
               print = TRUE
)
# hmm
# Find states with the Viterbi algorithm
v <- as.integer(factor(viterbi(hmm, counts), levels = hmm$StateNames))
plot(sqrt(counts), col = v)

# CBS
cna <- CNA(matrix(sqrt(counts), ncol = 1),
           chrom = as.character(seqnames(query)),
           maploc = start(query),
           data.type = "logratio",
           presorted = TRUE
)
scna <- segment(cna,
                undo.splits = "sdundo",
                undo.SD = 1.5,
                verbose = 1
)
seq <- with(scna$output, rep(seg.mean, num.mark))
# plot(scna)
# plot(seq)
q <- quantile(seq, .95)
seq2 <- pmin(seq, q)
# plot(seq2)
km <- kmeans(seq2, centers = 3)
query$states <- km$cluster
plot(sqrt(counts), col = km$cluster)

