options(stringsAsFactors = FALSE)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pander)
library(scatterplot3d)
library(limma)

setwd("~/Documents/Work/Teaching/BIOS691_Cancer_Bioinformatics/static/slides/09_Clustering/lab")
# Get the data from http://odin.mdacc.tmc.edu/~kdo/TeachBioinf/Projects%20&%20Data%20Sets/nci60.tsv
mtx <- read.delim("../data/nci60.tsv", sep = "\t", row.names = 1)
dim(mtx)
colnames(mtx)
boxplot(mtx)
# Create annotations 
groups <- sapply(colnames(mtx), function(x) { obj <- strsplit(x, ".", fixed = TRUE); unlist(obj)[[1]][1]}) %>% as.character()
groups[ groups == "K562A" ] <- "K562" # A couple of manual tweaks
groups[ groups == "K562B" ] <- "K562"
groups[ groups == "MCF7A" ] <- "MCF7" # A couple more
groups[ groups == "MCF7D" ] <- "MCF7"
groups %>% table

# PCA: Check for batch effects. Select one batch, to color points by its assignment
pca <- mtx %>% t %>% scale %>% prcomp 
pca <- prcomp(t(scale(mtx)))
data.frame(summary(pca)$importance)[, 1:5] %>% pander # Percent of variance explained

colorby <- "cells"
pt <- ggplot(data = data.frame(pca$x, cells = groups, samples = groups, stringsAsFactors = F), 
             aes(x = as.numeric(PC1), y = as.numeric(PC2), label = samples)) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  ggtitle(paste("PCA with batch, coloring by ", colorby)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
  geom_text_repel(colour = "black", size = 3) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  labs(color = colorby) +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca)$importance[2,2] * 100, digits = 2), "% variability" ))
plot(pt)
# ggsave(filename = "Figure.pdf", plot = pt, height = 8, width = 11)

# devtools::install_github("fawda123/ggord")
# library(ggord)
# ggord(pca, groups, arrow = NULL)

#install.packages("ggfortify")
library(ggfortify)
autoplot(prcomp(mtx %>% t %>% scale), data = data.frame(groups = groups), colour = "groups")

# library(MDmisc)
# library(pcaGoPromoter)
library(ellipse)
library(gridExtra)
library(grid)
library(ggplot2)
# pca_func(mtx, groups, title = "PCA")

# Let's introduce batch into a half of the matrix - raise it to the power of 2
mtx_with_batch <- cbind(mtx[, 1:round(ncol(mtx) / 2)], mtx[, (round(ncol(mtx) / 2) + 1):ncol(mtx)]^2)
boxplot(mtx_with_batch)
# Batch variable, known to us
batch <- c(rep(1, round(ncol(mtx) / 2)), rep(2, ncol(mtx) - round(ncol(mtx) / 2)))
batch

# PCA with batch
pca <- mtx_with_batch %>% scale %>% t %>% prcomp
data.frame(summary(pca)$importance)[, 1:5] %>% pander # Percent of variance explained

colorby <- "batch"
ggplot(data = data.frame(pca$x, cells = groups, samples = groups, batch = factor(batch), stringsAsFactors = F), 
       aes(x = as.numeric(PC1), y = as.numeric(PC2), label = samples)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
  geom_text(colour = "black", size = 3)

# A function to pull out p-value of LM. https://stackoverflow.com/questions/5587676/pull-out-p-values-and-r-squared-from-a-linear-regression
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
# Which covariate explains the PCs
covariates <- c("batch")
for (covariate in covariates){
  pca.lm <- lm( as.numeric(PC1) ~ factor(eval(parse(text = covariate))), data = data.frame(batch = batch, pca$x))
  print(paste(covariate, "accounts for", signif(summary(pca.lm)$adj.r.squared, 5), "variability explained by the first principle component, p-value", signif(lmp(pca.lm), 5)))
  pca.lm <- lm( as.numeric(PC2) ~ factor(eval(parse(text = covariate))), data = data.frame(batch = batch, pca$x))
  print(paste(covariate, "accounts for", signif(summary(pca.lm)$adj.r.squared, 5), "variability explained by the second principle component, p-value", signif(lmp(pca.lm), 5)))
  pca.lm <- lm( as.numeric(PC3) ~ factor(eval(parse(text = covariate))), data = data.frame(batch = batch, pca$x))
  print(paste(covariate, "accounts for", signif(summary(pca.lm)$adj.r.squared, 5), "variability explained by the third principle component, p-value", signif(lmp(pca.lm), 5)))
}

# How does it compare with MDS?
plotMDS(mtx_with_batch , col = ifelse(batch == 1, "red", "blue"), main = "Before batch removal")

# t-SNE
library(Rtsne)

tsne <- Rtsne(mtx %>% t %>% scale, dims = 2, perplexity=20, verbose=TRUE, max_iter = 500)

## Plotting
plot(tsne$Y, t='n', main="tsne")
colors <- rainbow(length(unique(groups)))
text(tsne$Y, labels = groups, col=colors[groups])

library(ggplot2)
library(ggrepel)
scores <- as.data.frame(tsne$Y)
rownames(scores) <- colnames(mtx)
colnames(scores) <- c("Comp.1", "Comp.2")
scores <- data.frame(scores, groups = groups)
pt <- ggplot(data=scores, aes(x=Comp.1, y=Comp.2, label=groups, color = groups)) +
  #  labs(x = "PC1 (23.06%)", y = "PC2 (7.95%)") +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  theme(legend.position = c(0.1, 0.12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(legend.key.size = unit(10, "mm")) +
  geom_point(aes(shape = groups), size = 3) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text_repel(colour = "black", size = 3)
plot(pt)

# PCD 3D
# http://davetang.org/muse/2015/02/12/animated-plots-using-r/
scores <- as.data.frame(pca$x)
scatterplot3d(scores[, c("PC1", "PC2", "PC3")], color=as.numeric(factor(batch)))

# The following would be MAC-specific
rename <- function(x){
  if (x < 10) {
    return(name <- paste('000',i,'plot.png',sep=''))
  }
  if (x < 100 && i >= 10) {
    return(name <- paste('00',i,'plot.png', sep=''))
  }
  if (x >= 100) {
    return(name <- paste('0', i,'plot.png', sep=''))
  }
}
frames <- 360
#loop through plots
for(i in 1:frames){
  name <- rename(i)
  #saves the plot as a .png file in the working directory
  png(name)
  s3d <- scatterplot3d(scores[, c("PC1", "PC2", "PC3")],
                       main=paste("Angle", i),
                       angle=i,
                       pch=19,
                       cex.symbols=0.5,
                       color=as.numeric(factor(batch)))
  s3d.coords <- s3d$xyz.convert(scores[, c("PC1", "PC2", "PC3")])
  text(s3d.coords$x, s3d.coords$y, labels=groups, pos=2, offset=0.5, cex=0.7)
  dev.off()
}
# The following is Mac/Unix-specific
# Conversion of multiple PNG files into one GIF using ImageMagik, https://www.imagemagick.org/script/index.php
my_command <- 'convert *.png -delay 1 -loop 0 3d.gif'
system(my_command)
system("rm *.png")
