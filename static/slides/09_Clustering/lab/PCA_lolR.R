# https://cran.r-project.org/web/packages/lolR/index.html
install.packages("lolR")
library(lolR)
library(ggplot2)
library(MASS)

n=400
d=30
r=3

testdat <- lol.sims.cigar(n, d)
X <- testdat$X
Y <- testdat$Y

data <- data.frame(x1=X[,1], x2=X[,2], y=Y)
data$y <- factor(data$y)
ggplot(data, aes(x=x1, y=x2, color=y)) +
  geom_point() +
  xlab("x1") +
  ylab("x2") +
  ggtitle("Simulated Data")

result <- lol.project.pca(X, r)

data <- data.frame(x1=result$Xr[,1], x2=result$Xr[,2], y=Y)
data$y <- factor(data$y)
ggplot(data, aes(x=x1, y=x2, color=y)) +
  geom_point() +
  xlab("x1") +
  ylab("x2") +
  ggtitle("Projected Data using PCA")

cor(result$Xr[,1], X)
cor(result$Xr[,2], X)

liney <- MASS::lda(result$Xr, Y)
result <- predict(liney, result$Xr)
lhat <- 1 - sum(result$class == Y)/length(Y)

data <- data.frame(x1=result$x[,1], y=Y)
data$y <- factor(data$y)
ggplot(data, aes(x=x1, fill=y)) +
  geom_density(adjust=1.5, alpha=0.6) +
  xlab("$x_1$") +
  ylab("Density") +
  ggtitle(sprintf("PCA-LDA, L = %.2f", lhat))
