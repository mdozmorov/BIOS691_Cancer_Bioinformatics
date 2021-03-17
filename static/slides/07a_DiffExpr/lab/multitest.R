### Two Sample T-test
group1 <- c(2013.7, 2141.9, 2040.2, 1973.3, 2162.2, 1994.8, 1913.3, 2068.7)
group2 <- c(1974.6, 2027.6, 1914.8, 1955.8, 1963, 2025.5, 1865.1, 1922.4)

length(group1)
length(group2)
mean(group1)
mean(group2)
var(group1)
var(group2)
t.value <- (mean(group1) - mean(group2))/sqrt(var(group1)/length(group1) + var(group2)/length(group2))
t.test(group1, group2)
t <- seq(-3.5, 3.5, by = 0.01)
?dt
plot(t, dt(t, df = 12.116))
pt(t.value, df = 12.116)
?pt
pt(t.value, df = 12.116, lower.tail = FALSE)
2 * pt(t.value, df = 12.116, lower.tail = FALSE)

## Alternatively, an equivalent way of performing the t-test is
length(group1)
length(group2)
group <- rep(0:1, each = 8)
all <- c(group1, group2)
t.test(all ~ group)
out <- t.test(all ~ group)
names(out)
out$p.value

### Multiple Hypothesis Testing Issues
data <- matrix(NA, ncol = 16, nrow = 10000)
for (i in 1:10000) {
  data[i, ] <- rnorm(16, 2000, 100)
}
data[1, ]
p.value <- numeric()
for (i in 1:10000) {
  p.value[i] <- t.test(data[i, ] ~ group)$p.value
}
hist(p.value)
sum(p.value < 0.05)

g <- 10000
alpha <- 0.05
bonferroni <- alpha/g
bonferroni
SidakSS <- 1 - (1 - alpha)^(1/g)
SidakSS
sum(p.value < bonferroni)
sum(p.value < SidakSS)
holm <- numeric()
SidakSD <- numeric()
for (i in 1:g) {
  holm[i] <- alpha/(g - i + 1)
  SidakSD[i] <- 1 - (1 - alpha)^(1/(g - i + 1))
}
holm
SidakSD
sum(sort(p.value) < holm)
sum(sort(p.value) < SidakSD)

# Real world example
library(Biobase)
load("/Users/mdozmorov/Documents/Work/Teaching/BIOS567.2017/assets/07a_DiffExpr/lab/Methyl.RData")
ls()
methyl # 1490 features, 95 samples
exprs(methyl)[1:5, 1:2]
pData(methyl)
plot(density(exprs(methyl)[, 1]))
### Now note there are 1490 CpG sites (features).

### FWER adjusted alpha levels
g <- 1490
alpha <- 0.05
bonferroni <- alpha/g
SidakSS <- 1 - (1 - alpha)^(1/g)
holm <- numeric()
SidakSD <- numeric()
for (i in 1:g) {
  holm[i] <- alpha/(g - i + 1)
  SidakSD[i] <- 1 - (1 - alpha)^(1/(g - i + 1))
}
cbind(holm, SidakSD)[1:10, ]

### For each CpG site, apply a two-sample t-test to test for differential
### expression between males and females
table(pData(methyl)$Gender.rep)
gender <- pData(methyl)$Gender.rep - 1
table(gender)
pvalue <- numeric()
for (i in 1:dim(exprs(methyl))[1]) {
  exprs <- as.numeric(exprs(methyl)[i, ])
  pvalue[i] <- t.test(exprs[gender == 0], exprs[gender == 1])$p.value
}
sum(pvalue < bonferroni)
sum(pvalue < SidakSS)
sum(sort(pvalue) < holm)
sum(sort(pvalue) < SidakSD)

### Benjamini & Hochberg method
k <- 1:g
sorted.p <- sort(pvalue)
bh <- ifelse(sorted.p <= (0.05 * k)/1490, TRUE, FALSE)
sum(bh)

#### p.adjust function
library(stats)
by <- p.adjust(pvalue, method = "BH")
sum(by < 0.05)

### Get the FDR using the default q-value method
library(qvalue)
qval <- qvalue(pvalue)
names(qval)
sum(qval$qvalues < 0.05)

############################################### FDR methods are available in the limma package but the results differ
############################################### slightly due to the moderated t-test
library(limma)
f <- factor(gender, levels = unique(gender))
design <- model.matrix(~0 + f)
# colnames(design)<-c('gender')
fit <- lmFit(methyl, design = design)
contrast.matrix <- makeContrasts(deg = f1 - f0, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit3 <- eBayes(fit2)
### adjust method can be either 'none', 'BH', 'BY' and 'holm'

############## none
out <- topTable(fit3, adjust.method = "none", number = Inf, sort.by = "none")
sign <- rownames(out)[out$adj.P.Val < 0.05]
length(sign)

############### BH
out <- topTable(fit3, adjust.method = "BH", number = Inf, sort.by = "none")
sign <- rownames(out)[out$adj.P.Val < 0.05]
length(sign)

############### Holm
out <- topTable(fit3, adjust.method = "holm", number = Inf, sort.by = "none")
sign <- rownames(out)[out$adj.P.Val < 0.05]
length(sign)

