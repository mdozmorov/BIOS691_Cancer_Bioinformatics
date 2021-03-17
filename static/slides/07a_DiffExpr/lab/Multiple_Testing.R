# This lab investigates the issues of multiple comparison associated with large P small N data (like microarray data)
# in particular the problems of dealing with correlated data
# The lab is a simulation of a microarray experiment and the analysis
# the simulated experment has two groups of N = 10  samples with M = 1,000 genes
# The value of doing the simulation is that you know the true values so you can
# compare the effectiveness of various procedures

M <- 1000 # number of 'genes' in the 'study'
N <- 10  # number of samples in each of two groups
gene.means <- rep(0,M) # for convenience set all gene means to be 0 
treatment.true.values <- control.true.values <- matrix( 0, gene.means, nr=M, nc=N)
   # true values are means 
rownames(treatment.true.values) <- rownames(control.true.values) <- paste("Gene",1:M)

# This section aims to challenge your intuitions about p-values 
# You'll generate some completely random data and find that, 
# even at a stringent p-value threshold, several genes appear 'significant'
noise1 <- matrix(rnorm(M*N*2) ,nr = M, nc=2*N) # generate random independent errors ('noise')
sim.data1 <- cbind( treatment.true.values, control.true.values) + noise1 # no real differences between two groups
library(genefilter) # we load this library for fastT to do many t-tests in parallel
t.scores <- fastT( sim.data1, 1:N, (N+1):(2*N))$z # generate many t-scores in parallel  
pvals <- 2 * pt( -abs(t.scores), df = 2*N-2) # compute p-values
# how do you think all these p-values should be distributed? why?
hist(t.scores)
hist(pvals)

# How many 'significant' changes are there?
sum( pvals < 0.01)
# how does the distribution of p-values look?
hist( pvals,breaks=20) 
# Try out the Bonferroni correction
sum( pvals < 0.05 / M ) # how many are even moderately significant after 'correction'?

# Now let's examine a situation in which we know there are some real effects
treatment.true.values[1:10,] <- 2 
# ten genes changed by 2 SD of 'noise' (normal variation within experimental groups)
# this is actually quite a large effect size; in a single test we'd have a 99% chance of detecting the difference with samples of 10 and 10
# there are no changes in the remaining 990 genes
sim.data1 <- cbind( treatment.true.values, control.true.values) + noise1 # now 10 real effects
t.scores <- fastT( sim.data1, 1:N, (N+1):(2*N))$z # generate many t-scores quickly  
# what's returned by fastT? try tmp <- fastT(sim.data1, 1:10, 11:20); names(tmp)
pvals <- 2 * pt( -abs(t.scores), df = 2*N-2)
# there should be more p-values under 0.01 ...
sum( pvals < 0.01) # we do get some apparently significant genes, out of 1,000...
sum( pvals[1:10] < 0.01) # ... but only some of them are the truly changed genes
sum( pvals[11:1000] < 0.01) # while others are false positives
# Try this with different thresholds 
for ( threshold in c(.001,.002,.005,.01,.02)) 
	cat('\ntotal significant:',sum( pvals < threshold),';\ttrue:', sum( pvals[1:10] < threshold),'\tfalse:', sum( pvals[11:1000] < threshold))

# this is a recurrent theme in multiple testing ... 
# we can't get all the changed genes without getting a flood of false positives
# These are unusually strong effects (2 SD); usually it's much worse; try effect size 1

# how does the distribution of p-values look now?
hist( pvals,breaks=50) # ... maybe a slight pile-up at the left end (p~0)
# Try out the Bonferroni correction
pvals.bonf <- pmin( M*pvals, 1 ) 
sum( pvals.bonf < .1 ) # how many are modestly significant?
sum( pvals.bonf[1:10] < .1 ) # how many 'false negatives' by Bonferroni?

## Optional: Try out Holm's and Hochberg's FWER procedure at alpha = 0.1
sum( sort(pvals) < 0.1/M:1 )  # probably won't find anything

# Now simulate many more data sets with the same number of truly changed genes
# but different random fluctuations and test the Bonferroni correction
K = 100 ; true.counts <- false.counts <- numeric(K)
treatment.true.values[1:10,] <- 2 # ten genes changed by 2 SD; no change in the remainder
for ( kk in 1:K ) {
  noise1 <- matrix(rnorm(M*N*2) ,nr = M, nc=2*N) # generate random independent errors ('noise')
  sim.data1 <- cbind( treatment.true.values, control.true.values) + noise1 # no real effects
  t.scores <- fastT( sim.data1, 1:N, (N+1):(2*N))$z 
  pvals <- 2 * pt( -abs(t.scores), df = 2*N-2)
  pvals.bonf <- pmin( M*pvals, 1 ) 
  true.counts[kk] <- sum( pvals.bonf[1:10] < .1)  
  false.counts[kk] <- sum(pvals.bonf[11:M] < .1)
  if ( kk < 10) cat( '\nBonferroni',kk,':', sum(pvals.bonf < .1),'; true', sum(pvals.bonf[1:10]<.1 ) )
}
# How often does the Bonferroni procedure give any false positives?
table( true.counts, false.counts)
# we never get all the truly changed genes, and we usually get a fair number of false positives

# We've been considering an idealized situation - the 'normal' variations of one gene are unrelated to all the others
# In reality genes are always somewhat co-regulated (and technical artifacts are correlated)
# So let's now try a few random data sets with correlated errors
for ( kk in 1:K ) {
  noise2 <- sqrt(1/2)*noise1 + sqrt(1/2) * matrix( rep(1,M),nr=M,nc=2*N) * matrix( rnorm(N*2),nr=M, nc=2*N,byrow=T) 
  	# the second term is a fluctuation common to all genes
#  noise2 <- sqrt(1/2)*noise1 + sqrt(1/2) * matrix( c(rep(1,M/2),rep(-1,M/2)),nr=M,nc=2*N) * matrix( rnorm(N*2),nr=M, nc=2*N,byrow=T) 
  	# the second term is a fluctuation opposite in sign for half the genes
  sim.data2 <- cbind( treatment.true.values, control.true.values) + noise2 # no real effects
  t.scores <- fastT( sim.data2, 1:N, (N+1):(2*N))$z 
  pvals <- 2 * pt( -abs(t.scores), df = 2*N-2)
  pvals.bonf <- pmin( M*pvals, 1 ) 
  true.counts[kk] <- sum( pvals.bonf[1:10] < .1)  
  false.counts[kk] <- sum(pvals.bonf[11:M] < .1)
#  cat( '\nBonferroni',kk,':', sum(pvals.bonf < .1),'; true', sum(pvals.bonf[1:10]<.1 ) )
}
# How do the numbers of positives and false positives differ from the IID case?
table( true.counts, false.counts)
plot(jitter( true.counts), jitter(false.counts)) # the 'jitter' spreads out points that would be over-plotted
# the number of false positives goes up with the number of true positives
# To look at negative correlations remove the comment from the second 'noise2' line above and run the loop again

# In real genomic data the correlation structure is much more complex than either of these examples
# The Bonferroni procedure is conservative and works for all these cases. We'll see that some other procedures aren't so reliable


## In real genomic data we are usually dealing with a range of real effect sizes
# (i.e. by how much the gene expressions are changed)
## Now let's simulate that more realistic situation: 
# set out a distribution of 'real' effects
mean.change <- c( 3 * (1 - 1:100/101), rep(0, M-100 ))
# effect sizes uniformly distributed from 0 to 3
# Plot to see what these effects are
plot(mean.change)
treatment.true.values[1:M,] <- mean.change

# Add random independent errors ('noise') to these true values
noise1 <- matrix( rnorm(M*N*2), nr = M, nc=2*N ) 
sim.data1 <- cbind( treatment.true.values, control.true.values) + noise1  

## Select genes by simple t-scores  ###
t.scores <- fastT( sim.data1, 1:N, (N+1):(2*N))$z # fastT does parallel t-tests in Fortran
pvals <- 2 * pt( -abs(t.scores), df = 2*N-2)
hist( pvals,breaks=50) # why does it look different from the one before?

# Estimate the expected proportion of false positives if you pick .01 as a threshold
0.01*M / sum( pvals < .01) # you'd expect 1% of all genes to have p < 0.01
# now compare the actual proportion of false positives at p < .01
sum( pvals[101:M] < .01) / sum( pvals < .01) 
# How close is the simple estimate?     
# Start again from 'noise1 <- ' and do this a few times to get a feel for how variable the proportion is

# OPTIONAL: Plot p-value .vs. real effect size to see the power of the test
plot( mean.change, -log(pvals), main='Power to detect change')
abline( h = -log( 1 / M ), col=2 ) # indicate the p-value needed to get about 1 false positive by Bonferroni
# Of the 'genes' with an effect size of 1, about what fraction of true positives do you miss?
# Plot estimated effect size against real effect size
estimated.mean.change <- apply( sim.data1[,1:10], 1, mean) - apply( sim.data1[,11:20], 1, mean)
plot( mean.change, estimated.mean.change ) ; abline( 0, 1, col=2)
# check that errors in estimated effect size are about what are predicted by the statistical theory ...
var ( mean.change - estimated.mean.change ) # What should this be? 

### Prelude to FDR
# Examine the trade-off between false positives to false negatives at various p-value thresholds for a simple data set and for a complex data set
counts.fp <- counts.fn <- counts.tp <- numeric()
for ( ii in 1:100 ) {  # for each ii set a threshold of p < exp( -ii/10)
  counts.fp[ii] <- sum( log(pvals[101:M]) < - ii / 10 )
  counts.fn[ii] <- sum( log(pvals[1:100]) > - ii/10 )
  counts.tp[ii] <- sum( log(pvals[1:100]) < - ii/10 )
}
# Construct an ROC (Receiver operating characteristic) plot
# This records how many false positives you need to accept in order to get a given fraction of true positives
plot( counts.fp / M, 1 - counts.fn/100, ylim=c(0,1), ylab='Sensitivity', xlab='False Positives',type='b'); abline(0,1,col=3)
points( counts.fp[10*1:10] / M, 1 - counts.fn[10*1:10]/100, pch=16,col=2 )
   # the red dots just help you count: every tenth threshold is red
# What p-value threshold seems to be the sweet spot?
# Now plot the proportion of false positives as a function of p-value threshold
total.counts <- counts.fp + counts.tp
plot( exp(-(1:100)/10),counts.fp/total.counts, ylim=c(0,1), log='x', ylab='False Positive Proportion',xlab='p-value threshold') 
# What p-value thresholds seem to give reasonable proportions of false positives?
# We'd like some way to pick a threshold on real data when we don't know the truth

############  FDR Procedures   ################
M <- 1000 # number of 'genes' in the 'study'
N <- 10  # number of samples in each of two groups
gene.means <- rep(0,M) # for convenience set all gene means to be 0 
treatment.true.values <- control.true.values <- matrix( 0, gene.means, nr=M, nc=N)
rownames(treatment.true.values) <- rownames(control.true.values) <- paste("Gene",1:M)
treatment.true.values[1:10,] <- 2 # ten genes changed by 2 SD; no change in the remainder
noise1 <- matrix(rnorm(M*N*2) ,nr = M, nc=2*N) # generate random independent errors ('noise')
sim.data1 <- cbind( treatment.true.values, control.true.values) + noise1 # no real differences between two groups
library(genefilter) # we load this library for fastT to do many t-tests in parallel
t.scores <- fastT( sim.data1, 1:N, (N+1):(2*N)) # generate many t-scores in parallel  
pvals <- 2 * pt( -abs(t.scores$z), df = 2*N-2) # compute p-values
# Try out the Benjamini-Hochberg FDR procedure
alpha = 0.1   # Aim for 10% false discoveries
# how many genes would be dclared significant?
satisfy.set <- which( sort(pvals) < alpha * 1:M / M ) # the set of p-values that satisfies the condition 
(k<-max(satisfy.set)) ; sort(pvals)[k]
# select the 'genes' with p-values up to and including largest p-value satisfying criterion
selected <- order( pvals )[1:k] # 'genes' with smallest to k-th p-values
# How many false positives did you actually get?
sum( selected %in% 11:1000 ) / length(selected) 
# how well does this compare with your target proportion?

# try this again a few times (starting from noise1 <-) to see how variable the B-H procedure can be on IID data

# Construct Q-values
Q <- numeric(); Q[order(pvals)] <- M / 1:M * sort(pvals )
for ( ii in M:2 ) Q[ order(pvals)[ii-1] ] <- min( Q[ order(pvals)[(ii-1):ii] ] )
# Show the relationship between apparent change and Q-value in a 'Volcano Plot'
estimated.mean.change <- apply( sim.data1[,1:10], 1, mean) - apply( sim.data1[,11:20], 1, mean)
plot( estimated.mean.change, -log(Q))
                                                    
# Now examine how the same procedures work with correlated errors and variable effect sizes
# as occurs in real life
# Construct correlated errors where one factor accounts for 50% of variance across all genes
# NB this is a little bit more than usual 
mean.change <- c( 3 * (1 - 1:100/101), rep(0, M-100 ))
treatment.true.values[1:M,] <- mean.change
noise1 <- matrix(rnorm(M*N*2) ,nr = M, nc=2*N) 
noise2 <- sqrt(1/2)*noise1 + sqrt(1/2) * matrix( rep(1,M),nr=M,nc=2*N) * matrix( rnorm(N*2),nr=M, nc=2*N,byrow=T) 
#noise2 <- sqrt(1/2)*noise1 + sqrt(1/2) * matrix( c(rep(1,M/2),rep(-1,M/2)),nr=M,nc=2*N) * matrix( rnorm(N*2),nr=M, nc=2*N,byrow=T) 
var(c(noise1)); var(c(noise2)) # check that the errors are comparable in size to previous situation
sim.data2 <- cbind( treatment.true.values, control.true.values) + noise2  
t.scores <- fastT( sim.data2, 1:N, (N+1):(2*N))$z
pvals <- 2 * pt( -abs(t.scores), df = 2*N-2 )
plot( mean.change, -log10(pvals))
hist(pvals,breaks=50) # notice that there appear to be 'waves' in the p-value histogram
satisfy.set <- which( sort(pvals) < alpha * 1:M / M ) # the set of p-values that satisfies the condition 
(k<-max(satisfy.set)) ; sort(pvals)[k]
# select the 'genes' with p-values up to and including largest p-value satisfying criterion
selected <- order( pvals )[1:k] # 'genes' with smallest to k-th p-values
# How many false positives did you actually get?
sum( selected %in% 11:1000 ) / length(selected) 
# How does this variation compare with the situation with independent errors?



# investigate the variability in the numbers of false positives systematically under different correlation structures
# you'll have to uncomment the appropriate lines below
counts.p.Q <- counts.fp.Q <- Q <- numeric()  # numbers of p-values under 0.01 and false positives
for ( ii in 1:100 ) {
  noise1 <- matrix(rnorm(M*N*2) ,nr = M, nc=2*N) 
  # positively correlated errors
  noise2 <- sqrt(2/3)*noise1 + sqrt(1/3)* matrix( c(rep(1,M/2),rep(0,M/2)),nr=M,nc=2*N) * matrix( rnorm(N*2),nr=M, nc=2*N,byrow=T) 
  # negatively correlated errors
#  noise3 <- sqrt(2/3)*noise1 + sqrt(1/3)* matrix( c(rep(1,M/2),rep(-1,M/2)),nr=M,nc=2*N) * matrix( rnorm(N*2),nr=M, nc=2*N,byrow=T) 
  sim.data2 <- cbind( treatment.true.values, control.true.values) + noise2  
  t.scores <- fastT( sim.data2, 1:N, (N+1):(2*N))$z
  pvals <- 2 * pt( -abs(t.scores), df = 2*N-2 )
  Q[order(pvals)] <- cummin( (M/1:M * sort( pvals))[M:1])[M:1]    # Q values for FDR
  counts.p.Q[ii] <- sum( Q < .1 )
  counts.fp.Q[ii] <- sum( Q[101:M] < .1) 
}

# Compare the simple estimated FDR to the actual FDR
mean( counts.fp.Q / counts.p.Q )  # what should this be?
hist( counts.fp.Q / counts.p.Q )
#sd( counts.fp.Q / counts.p.Q ) / mean( counts.fp.Q / counts.p.Q ) # coefficient of variation
# What are the real FDR's for selection at a nominal 10% Q value?
plot( counts.p.Q , counts.fp.Q / counts.p.Q, xlab = 'Number significant', ylab='Actual FDR' ) ; abline(h=0.1,col=2); abline(h=0)
# do these results surprise you?
# How confident can you be in the FDR estimates from B-H

# Storey Procedure
# work with the same data set
M <- 1000 # number of 'genes' in the 'study'
N <- 10  # number of samples in each of two groups
treatment.true.values <- control.true.values <- matrix( 0, nr=M, nc=N)
   # true values are means 
rownames(treatment.true.values) <- rownames(control.true.values) <- paste("Gene",1:M)
mean.change <- c( 3 * (1 - 1:100/101), rep(0, M-100 ))
treatment.true.values[1:M,] <- mean.change
noise1 <- matrix(rnorm(M*N*2) ,nr = M, nc=2*N) 
  # positively correlated errors
noise2 <- sqrt(2/3)*noise1 + sqrt(1/3)* matrix( c(rep(1,M/2),rep(0,M/2)),nr=M,nc=2*N) * matrix( rnorm(N*2),nr=M, nc=2*N,byrow=T) 
sim.data2 <- cbind( treatment.true.values, control.true.values) + noise2  
t.scores <- fastT( sim.data2, 1:N, (N+1):(2*N))$z
pvals <- 2 * pt( -abs(t.scores), df = 2*N-2 )

pFDR <- function( pvals, p1 ) {
	# Step 1: estimate pi0
	pi0.est <- 2*mean( pvals > .5)
	M <- length(pvals) 
	# Step 2 adjust naive estimate
	pi0.est * p1 * M / max( sum(pvals<p1), 1 ) /(1-(1-p1)^M) 
}
K=1000 # bootstrap samples
bs<-numeric(K)
for (kk in 1:K) {
	bs[kk]<-pFDR( pvals[sample(1:length(pvals),replace=T)], .01)
}
# HW assignment ... run simulations like those in the demo to compare pFDR and B-H FDR estimates and to test how accurate they are under conditions of i) independence, ii) positive and iii) negative correlation. Where are they most different? (Bonus) Devise your own correlation schemes. Is the bootstrap confidence interval for pFDR accurate?

# Illustrate limma
source("http://bioconductor.org/biocLite.R")
biocLite("limma")
library(limma)
limmaUsersGuide(view=TRUE) # ?limma doesn't help much
design <- model.matrix(~ 0+factor(c(rep(1,10),rep(2,10))))
design 
colnames(design) <- c("group1", "group2")
fit <- lmFit(sim.data2, design)
summary(fit)
mean(fit$sigma^2)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
plot( fit$sigma^2, fit2$s2.post,xlim=c(0,2),ylim=c(0,2));abline(0,1,col=2)
topTable(fit2, coef=1, adjust="BH")
plot(t.scores,fit2$t); abline(0,-1,col=2)        
plot( mean.change, -t.scores)
points(mean.change+0.05, fit2$t, col=2)
