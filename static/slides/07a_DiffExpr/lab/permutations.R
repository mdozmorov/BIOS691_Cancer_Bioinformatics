# Total genes
a <- c(rep("On Other Chromosomes", 160), rep("On Chromosome 1", 40)) 
total.all   <- length(a)                        # Total genes
total.chr   <- sum(a == "On Chromosome 1")      # genes on a specific chromosome
total.other <- sum(a == "On Other Chromosomes") # genes on other chromosomes
all.equal(total.all, total.chr + total.other)   # they sum up to total

b <- c(rep("On Other Chromosomes", 35), rep("On Chromosome 1", 15)) # Selected genes (DEGs)
selected.chr   <- sum(b == "On Chromosome 1")
selected.other <- sum(b == "On Other Chromosomes")

# Obsserved 2x2 contingency table
mtx <- matrix(c(selected.chr, selected.other, total.chr - selected.chr, total.other - selected.other), 
              ncol = 2, dimnames = list(c("On Chromosome 1", "On Other Chromosomes"), 
                                        c("DEGs", "not DEGs")))
mtx                      # See how the table looks like
sum(mtx)                 # Make sure the total equals the total number of genes
fisher.test(mtx)$p.value # Fisher's exact p-value
chisq.test(mtx)$p.value  # Chi-square exact p-value

# Permuatation p-values
n <- 1:10000 # Number of permutations

# Permutation p-value using test statistics
estimate <- chisq.test(mtx)$statistic # Observed test statistics

set.seed(2)
# A vector of permutation t-statistics
estimate.rand <- sapply(n, function(x) {
  b.rand <- sample(a, size = length(b))                        # Randomly select genes
  selected.chr.rand <- sum(b.rand == "On Chromosome 1")        # How many are on chromosome
  selected.other.rand <- sum(b.rand == "On Other Chromosomes") # How many are not on chromosome
  # Permutation 2x2 contingency table
  mtx.rand <- matrix(c(selected.chr.rand, selected.other.rand, total.chr - selected.chr.rand, total.other - selected.other.rand), ncol = 2)
  chisq.test(mtx.rand)$statistic # Permutation test statistics
})
# Permutation p-value
sum(estimate.rand >= estimate) / length(n)

# Permutation p-value using difference in counts. H0: expected - observed = 0
estimate <- selected.chr # Observed counts

set.seed(1)
n <- 1:10000
# A vector of expected counts
estimate.rand <- sapply(n, function(x) {
  b.rand <- sample(a, size = length(b))
  sum(b.rand == "On Chromosome 1")
})

sum((estimate.rand - estimate) >= 0) / length(n)
