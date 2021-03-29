# kmeans clustering continued
weight.kg <- c(15, 49, 13, 45, 85, 66, 12, 10)
height.cm <- c(95, 156, 95, 160, 178, 176, 90, 78)
data <- matrix(cbind(weight.kg, height.cm), nrow = 8)
dimnames(data)[[2]] <- c("weight.kg", "height.cm")

# K-means using built-in capabilities
?kmeans
kmeans.clust <- kmeans(data, centers = 2)
kmeans.clust
kmeans.clust$centers
kmeans.clust$cluster
plot(weight.kg, height.cm, type = "n")
text(weight.kg, height.cm, labels = kmeans.clust$cluster)
points(kmeans.clust$centers, pch = c(24, 19))

### PAM - partitioning around medoids
# Select samples to be centroids. The plot is for illustration only
plot(weight.kg, height.cm, type = "n")
text(weight.kg[c(3, 1, 2, 8)], height.cm[c(3, 1, 2, 8)], labels = 1)
text(weight.kg[c(4, 5, 6, 7)], height.cm[c(4, 5, 6, 7)], labels = 2)

library(cluster)
?pam
pam.out <- pam(data, k = 2)
names(pam.out)
pam.out$medoids
pam.out$clustering
summary(pam.out)
plot(data)
points(pam.out$medoids, pch = 10)
par(mfrow = c(1, 2))
plot(pam.out)
?plot.partition

### Identifying k
x1 <- c(rnorm(25, 0, 1), rnorm(25, 3, 1))
x2 <- c(rnorm(25, 3, 1), rnorm(25, 0, 1))
data <- cbind(x1, x2)
par(mfrow = c(1, 1))
plot(data)
pam.2 <- pam(data, k = 2)
pam.3 <- pam(data, k = 3)
pam.4 <- pam(data, k = 4)
pam.5 <- pam(data, k = 5)
names(pam.2)
pam.2$silinfo
pam.2$silinfo$avg.width
plot(2:5, c(pam.2$silinfo$avg.width, pam.3$silinfo$avg.width, pam.4$silinfo$avg.width, pam.5$silinfo$avg.width))


### Filtering and clustering the most variable genes/samples
dir(pattern = "RData")
load("BreastCancer.RData")
ls()
Breast.MAS5

### For an Affymetrix dataset, first filter out the control probe sets
controls <- grep("AFFX", featureNames(Breast.MAS5))
length(controls)
featureNames(Breast.MAS5)[controls]
Breast.MAS5 <- Breast.MAS5[-controls, ]  ### Exclude all Affy control probe sets
Breast.MAS5  ### reduced by 68

coef.var <- function(x) {
  sqrt(var(x))/mean(x)  ### Coefficient of variation
}

cv <- apply(exprs(Breast.MAS5), 1, coef.var)  # Variation across each row
length(cv)
summary(cv)
quantile(cv, 0.99)
filtered.breast <- Breast.MAS5[cv > quantile(cv, 0.99), ] # Keep top variable
filtered.breast

hclust <- hclust(dist(exprs(filtered.breast)), method = "average")
plot(hclust)

hclust <- hclust(dist(t(exprs(filtered.breast))), method = "average")
plot(hclust, cex = 0.5)

# Clustering based on correlation
hclust.rho <- hclust(as.dist((1 - cor(exprs(filtered.breast))) / 2), method = "average")
plot(hclust.rho, cex = 0.5)

hc1 <- hclust(as.dist(1 - abs(cor(exprs(filtered.breast)))), method = "average")
plot(hc1, cex = 0.5)

# Clustering correlations among genes
hc2 <- hclust(as.dist((1 - cor(t(exprs(filtered.breast)))) / 2), method = "average")
plot(hc2, cex = 0.5)

dend1 <- as.dendrogram(hc1) # Dendrogram object of clustered samples, contains order of clustered samples
dend2 <- as.dendrogram(hc2) # Dendrogram object of clustered genes, contains order of clustered genes
# to get red-green heatmap images
# source('http://www.bioconductor.org/biocLite.R')
# biocLite('gplots')
library(gplots)
heatmap(exprs(filtered.breast), col = greenred(75), Rowv = dend2, Colv = dend1, xlab = "Sample")
