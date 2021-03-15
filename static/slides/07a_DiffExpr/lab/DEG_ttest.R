rm(list = ls())
# read in CSV formatted data
# e=read.csv("http://rafalab.jhsph.edu/688/expression.csv",row.names=1)
setwd("/Users/mdozmorov/Documents/Work/Teaching/BIOS691_Cancer_Bioinformatics/static/slides/07a_DiffExpr/lab")
e <- read.csv("data/expression.csv.gz", row.names = 1)
e=log2(e) # log2 transform

m1=rowMeans(e[,1:3]) # group 1 
m2=rowMeans(e[,4:6]) # group 2
a=rowMeans(e)


which.min(m2-m1)
grp = rep(1:2, each=3)

o=order(m2-m1)
head(cbind(m1[o], m2[o]))

### installing packages
source("http://bioconductor.org/biocLite.R")
biocLite("genefilter")
library(genefilter)


s1=rowSds(e[,1:3])
s2=rowSds(e[,4:6])
ttest=(m2-m1)/sqrt(s1^2/3 + s2^2/3)

pval=2*(1-pt(abs(ttest),4))

plot(m2-m1, -log10(pval), main="volcano",xlab="Effect size",ylab="-log10(pvalue)")
abline(v = 2)
abline(v = -2)
abline(h = 2, col = "red")

smoothScatter(m2-m1, -log10(pval), main="volcano",xlab="Effect size",ylab="-log10(pvalue)")

tt = rowttests(as.matrix(e),factor(grp))

plot(density(m1))
lines(density(m2), col = "red")
legend("topleft", c("grp1","grp2"), col=1:2, pch = 15)

plot(density(e[,1]), ylim = c(0,0.2))
for(i in 2:6) lines(density(e[,i]), col = i)
legend("topleft", legend=1:6, col = 1:6, pch = 15)

