# Example taken from the nullranges package
# https://github.com/nullranges/nullranges
# BiocManager::install("nullranges/nullranges")
library(nullranges)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)

edb <- EnsDb.Hsapiens.v86
g <- genes(edb, filter=AnnotationFilterList(GeneIdFilter("ENSG", "startsWith")))
x <- keepSeqlevels(g,as.character(seq_len(12)),pruning.mode="coarse")

Ls=1e6
n=3
seg<-segment_density(g1,n=3,Ls=Ls,type = "HMM",boxplot = TRUE)

query <- tileGenome(seqlengths(x)[seqnames(x)@values], tilewidth = Ls, cut.last.tile.in.chrom = TRUE)
counts <- countOverlaps(query, x)
eps <- rnorm(length(counts), 0, .2)
plot(sqrt(counts) + eps)

# CBS
library(DNAcopy)
cna <- CNA(matrix(sqrt(counts) + eps, ncol = 1),
           chrom = as.character(seqnames(query)), # wont work for X,Y,MT
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
km <- kmeans(seq2, n)
query$states <- km$cluster
plot(sqrt(counts) + eps, col = km$cluster)

# HMM
library(RcppHMM)
hmm <- initPHMM(n)
hmm <- learnEM(hmm,
               counts,
               iter = 400,
               delta = 1e-5,
               print = TRUE
)
hmm
v <- as.integer(factor(viterbi(hmm, counts), levels = hmm$StateNames))
plot(sqrt(counts) + eps, col = v)






