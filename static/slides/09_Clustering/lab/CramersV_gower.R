library(dplyr)
set.seed(1)
# Matrix with four categories in each column
mtx <- data.frame(a=factor(round(runif(200, 0, 3))),
                  b=factor(round(runif(200, 0, 3))),
                  c=factor(round(runif(200, 0, 3))))

# To test if the order of categories matter, repeat clustering with releveled matrix.
# Relative order of categories should be the same in each column, e.g. category "2" must mean same thing everywhere
mtx <- lapply(mtx, function(x) relevel(x, "3")) %>% as.data.frame


# Function for calculating Cramer's V manually
cramer <- function(x, y){
  K <- nlevels(x)
  L <- nlevels(y)
  n <- length(x)
  chi2 <- chisq.test(x, y, correct = FALSE)
  # print(chi2$statistic)
  v <- as.numeric(sqrt(chi2$statistic/(n * min(K - 1, L - 1))))
  return(v)
}

# Compare manual and builtin functions
library(lsr)
library(profvis)

profvis({
  m1 <- cramer(mtx$a, mtx$b)
  
  m2 <- cramersV(mtx$a, mtx$b)
})

# Clustering, repeat with releveled matrix
sim <- matrix(1, nrow = ncol(mtx), ncol = ncol(mtx)) # Similarity matrix
rownames(sim) <- colnames(mtx)
colnames(sim) <- colnames(mtx)
for (i in 1:(nrow(sim) - 1)) {
  for (j in (i + 1):ncol(sim)) {
    y <- mtx[, i]
    x <- mtx[, j]
    sim[i, j] <- cramersV(x, y)
    sim[j, i] <- sim[i, j]
  }
}
# Distance matrix
dissim <- as.dist(1 - sim)
tree <- hclust(dissim, method = "ward.D")
plot(tree) 

# Trying Gower
library(cluster)
library(gplots)
library(dplyr)

# Data as factors
set.seed(1)
mtx_merged <- data.frame(a=as.factor(round(runif(200, 0, 3))),
                         b=as.factor(round(runif(200, 0, 3))),
                         c=as.factor(round(runif(200, 0, 3))), stringsAsFactors = F) %>% t %>% as.data.frame
dissim <- daisy(mtx_merged, metric="gower")  
as.matrix(dissim)
tree <- hclust(dissim, method = "ward.D")
heatmap.2(as.matrix(dissim), Rowv=as.dendrogram(tree), Colv=as.dendrogram(tree), trace="none", scale="none")

# Data as characters
set.seed(1)
mtx_merged <- data.frame(a=as.character(round(runif(200, 0, 3))),
                         b=as.character(round(runif(200, 0, 3))),
                         c=as.character(round(runif(200, 0, 3))), stringsAsFactors = F) %>% t %>% as.data.frame
dissim <- daisy(mtx_merged, metric="gower")  
as.matrix(dissim)
tree <- hclust(dissim, method = "ward.D")
heatmap.2(as.matrix(dissim), Rowv=as.dendrogram(tree), Colv=as.dendrogram(tree), trace="none", scale="none")

# A function to manually calculate gower distance
gower.dissim <- function(df){
  sim <- matrix(1, nrow = ncol(df), ncol = ncol(df)) # Similarity matrix
  rownames(sim) <- names(df)
  colnames(sim) <- names(df)
  
  for (i in 1:nrow(sim)) {
    for (j in 1:ncol(sim)) {
      d_ij <- (sum(df[, j] == df[, i], na.rm=TRUE)) / sum(!(is.na(df[,i])|is.na(df[,j])))
      sim[i, j] <- sim[j, i] <- d_ij
    }
  }
  # Distance matrix
  return(1 - sim)
} ##end function

dissim <- gower.dissim(t(mtx_merged))
tree <- hclust(as.dist(dissim), method = "ward.D")
heatmap.2(as.matrix(dissim), Rowv=as.dendrogram(tree), Colv=as.dendrogram(tree), trace="none", scale="none")
