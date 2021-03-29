#useful functions for psychometrics and personality research
#a growing compilation of simple functions that I have found useful for doing personality research
#last updated Feb  14, 2005
#by William Revelle
#includes 
#alpha.scale     #find coefficient alpha for a scale and a dataframe of items
#describe        give means, sd, skew, n, and se 
#summ.stats      #basic summary statistics by a grouping variable
#error.crosses   (error bars in two space)
#skew            find skew
#panel.cor       taken from the examples for pairs
#pairs.panels    adapted from panel.cor  --   gives a splom, histogram, and correlation matrix
#multi.hist  #plot multiple histograms
#correct.cor    #given a correlation matrix and a vector of reliabilities, correct for reliability
#fisherz        #convert pearson r to fisher z
#paired.r       #test for difference of dependent correlations
#count.pairwise  #count the number of good cases when doing pairwise analysis
#eigen.loadings  #convert eigen vector vectors to factor loadings by unnormalizing them
#principal       #yet another way to do a principal components analysis -- brute force eignvalue decomp 
#factor.congruence #find the factor congruence coeffiecints
#factor.model    #given a factor model, find the correlation matrix
#factor.residuals #how well does it fit?
#factor.rotate    # rotate two columns of a factor matrix by theta (in degrees)
#phi2poly       #convert a matrix of phi coefficients to polychoric correlations

#define a function to calculate coefficient alpha
alpha.scale=function (x,y)   #find coefficient alpha given a scale and a data.frame of the items in the scale
	{
		n=length(y)          #number of variables
		Vi=sum(diag(var(y,na.rm=TRUE)))     #sum of item variance
		Vt=var(x,na.rm=TRUE)                #total test variance
		((Vt-Vi)/Vt)*(n/(n-1))}              #alpha


#general descriptive statistics 
describe  <-  function (x, digits = 2,na.rm=TRUE)   #basic stats after dropping non-numeric data
{                                   #first, define a local function
    valid <- function(x) {      
        return(sum(!is.na(x)))
   		 }
    if (is.vector(x) )          #do it for vectors or 
    	{
    	    stats = matrix(rep(NA,6),ncol=6)    #create a temporary array
    
			stats[1, 1] <-  mean(x, na.rm=na.rm )
			stats[1, 2] <-  median(x,na.rm=na.rm  )
			stats[1, 3] <-  min(x, na.rm=na.rm )
			stats[1, 4] <-  max(x, na.rm=na.rm )
			stats[1, 5] <-  skew(x,na.rm=na.rm  )
			stats[1, 6] <-  valid(x )
        	len <- 1;
    	}
    	
    	
    else  {len = dim(x)[2]     #do it for matrices or data.frames
    
    stats = matrix(rep(NA,len*6),ncol=6)    #create a temporary array
    for (i in 1:len) {
    	if (is.numeric(x[,i])) {   #just do this for numeric data
			stats[i, 1] <-  mean(x[,i], na.rm=na.rm )
			stats[i, 2] <-  median(x[,i],na.rm=na.rm  )
			stats[i, 3] <-  min(x[,i], na.rm=na.rm )
			stats[i, 4] <-  max(x[,i], na.rm=na.rm )
			stats[i, 5] <-  skew(x[,i],na.rm=na.rm  )
			stats[i, 6] <-  valid(x[,i] )
        		}
    	}
    	}
   temp <-  data.frame(n = stats[,6],mean = stats[,1], sd = sd(x,na.rm=TRUE), median = stats[, 
        2],min= stats[,3],max=stats[,4], range=stats[,4]-stats[,3],skew = stats[, 5])
    answer <-  round(data.frame(V=seq(1:len),temp, se = temp$sd/sqrt(temp$n)),  digits)
     return(answer)
}


# find basic summary statistics by groups      #adapted from "Kickstarting R"
 summ.stats <- function (x,y) {               #data are x, grouping variable is y
  means <- tapply(x, y, mean, na.rm=TRUE)
  sds  <- tapply(x,y, sd, na.rm = TRUE)
  valid <- function (x) {return(sum(!is.na(x)) )}
  ns <- tapply(x,y, valid )
  se=sds/sqrt(ns)
  answer <-data.frame(mean=means,sd=sds,se=se,n=ns)
  }
  
#error bars in a two dimensional plot            
old.error.crosses <-function (x,y,z)  # x  and y are data frame with z points
    {for (i in 1:z)  
    	{xcen <- x$mean[i]
    	 ycen <- y$mean[i]
    	 xse  <- x$se[i]
    	 yse <-  y$se[i]
    	 arrows(xcen-xse,ycen,xcen+xse,ycen,length=.2, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
    	 arrows(xcen,ycen-yse,xcen,ycen+yse,length=.2, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
    	}	
   }

#error bars in a two dimensional plot with labels            
error.crosses <-function (x,y,labels=NULL,pos=NULL,arrow.len=.2,...)  # x  and y are data frame with 
    {z <- dim(x)[1]
     if (length(pos)==0) {locate<-rep(1,z)} else {locate <- pos}
     if (length(labels)==0) lab<- rep("",z) else lab <-labels
        for (i in 1:z)  
    	{xcen <- x$mean[i]
    	 ycen <- y$mean[i]
    	 xse  <- x$se[i]
    	 yse <-  y$se[i]
    	 arrows(xcen-xse,ycen,xcen+xse,ycen,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
    	 arrows(xcen,ycen-yse,xcen,ycen+yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
    	text(xcen,ycen,labels=lab[i],pos=positions[i],cex=1,offset=arrow.len+1)     #puts in labels for all points
    	}	
   }
#examples
# plot(mPA,mNA,pch=symb[condit],cex=4.5,col=colors[condit],bg=colors[condit],main="Postive vs. Negative Affect",xlim=c(-1,1.5),ylim=c(-1,1.5),xlab="Positive Affect",ylab="Negative Affect")
    
#error.crosses(paf.stats,naf.stats,4)     #put in x and y error bars! 


#find a vector of skews for a data matrix

 skew= function (x, na.rm = FALSE) 
 {
    if (na.rm)    x <- x[!is.na(x)]             #remove missing values
    sum((x - mean(x))^3)/(length(x) * sd(x)^3)  #calculate skew   
 }


     
 #some fancy graphics   -- adapted from help.cor

 panel.cor.scale <- function(x, y, digits=2, prefix="", cex.cor)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r = (cor(x, y,use="pairwise"))
         txt <- format(c(r, 0.123456789), digits=digits)[1]
         txt <- paste(prefix, txt, sep="")
         if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex * abs(r))
     }
     
  panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r = (cor(x, y,use="pairwise"))
         txt <- format(c(round(r,digits), 0.123456789), digits=digits)[1]
         txt <- paste(prefix, txt, sep="")
         if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex )
     }
     
     panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

pairs.panels <- function (x,y,smooth=TRUE,scale=FALSE,digits=2) #combines a splom, histograms, and correlations
      {if (smooth ){
         if (scale) {
             pairs(x,diag.panel=panel.hist,upper.panel=panel.cor.scale,lower.panel=panel.smooth)
                    }
                    else {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.smooth)
                    } #else  {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.smooth)
                    }
                   
                    else      #smooth is not true
             { if (scale) {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor.scale)
               } else  {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor) }
            } #end of else (smooth)
         
      }   #end of function
      
      
     
     
 multi.hist <- function(x,...) {nvar <- dim(x)[2]  #number of variables
     nsize=trunc(sqrt(nvar))+1   #size of graphic
     old.par <- par(no.readonly = TRUE) # all par settings which can be changed
     par(mfrow=c(nsize,nsize))       #set new graphic parameters
     for (i in 1:nvar) {
     name=names(x)[i]                #get the names for the variables
     hist(x[,i],main=name,xlab=name,...) }  #draw the histograms for each variable
     on.exit(par(old.par))   #set the graphic parameters back to the original
     }
     
   
   #function to do lowess fitting if we have missing data
  lowess.na <- function(x, y = NULL, f = 2/3,...) {  #do lowess with missing data
  
  	x1 <- subset(x,(!is.na(x)) &(!is.na(y)))
  	y1 <- subset(y, (!is.na(x)) &(!is.na(y)))
  lowess.na <- lowess(x1,y1,f, ...)
  }    
 
     
     
     
#function to replace upper triangle of matrix with unattenuated correlations, and the diagonal with reliabilities
#input is a correlation matrix and a vector of reliabilities
     
 correct.cor <- function(x,y) { n=dim(x)[1]   
        { x[1,1] <- y[1,1]
        for (i in 2:n) {
           x[i,i] <- y[1,i]   #put reliabilities on the diagonal
           k=i-1
           for (j in 1:k) {
              x[j,i] <- x[j,i]/sqrt(y[1,i]*y[1,j])  }   #fix the upper triangular part of the matrix
             }
           return(x)  }}
          
          
     
 #difference of dependent (paired) correlations (following Steiger,J., 1980, Tests for comparing elements of a correlation matrix. Psychological Bulletin, 87, 245-251)
 paired.r <- function(xy,xz,yz,n) {
       diff <- xy-xz
       determin=1-xy*xy - xz*xz - yz*yz + 2*xy*xz*yz
       av=(xy+xz)/2
       cube= (1-yz)*(1-yz)*(1-yz)
       t2 = diff * sqrt((n-1)*(1+yz)/(((2*(n-1)/(n-3))*determin+av*av*cube)))
       return(t2)
        }
      

  fisherz <- function(rho)  {0.5*log((1+rho)/(1-rho)) }   #converts r to z  
  
  
  #a function to call John Fox's polychor multiple times to form a matrix of polychoric correlations   
#needs package polychor

polychor.matrix <-function(x,y=NULL) {

 sizex <- dim(x)[2]
if (((is.data.frame(y))|(is.matrix(y))))  sizey<-dim(y)[2] 
   else  sizey <- dim(x)[2]
        
 result<-matrix(1,nrow=sizey,ncol=sizex)   #create the output array

 xnames<- names(x)
 colnames(result)<- names(x)

if (((is.data.frame(y))|(is.matrix(y))))  rownames(result) <- names(y)
    else  rownames(result) <- names(x)
 
    
if (!((is.data.frame(y))|(is.matrix(y)))) {     #default case returns a square matrix
   for (i in 2: sizex ) {
     for (j in 1:( i-1)) {
      result[j,i]<-polychor(table(x[,j],x[,i]) )
       result[i,j] <- result[j,i]
       }
      }
      }
     else {                   #if y is input, then return the rectangular array
        for (i in 1: sizex ) {
            for (j in 1:sizey) {
                 result[j,i]<-polychor(table(x[,i],y[,j]) )
                     }
                } }
   return (result) }
       
#count.pairwise   -- count the number of good cases when doing pairwise analysis (e.g.,for correlations)       
count.pairwise <- function(x,y) {
   sizex <- dim(x)[2]
if (is.data.frame(y))  sizey<-dim(y)[2] 
   else  sizey <- dim(x)[2]
        
 result<-matrix(1,nrow=sizey,ncol=sizex)   #create the output array

 xnames<- names(x)
 colnames(result)<- names(x)
 
 if (((is.data.frame(y))|(is.matrix(y))))  rownames(result) <- names(y)
    else  rownames(result) <- names(x)
 
    
if (!is.data.frame(y)) {     #default case returns a square matrix
   for (i in 2: sizex ) {
     for (j in 1:( i-1)) {
     result[j,i]<- sum((!is.na(x[,j]))&(!is.na(x[,i])))
     result[i,j] <-  result[j,i]
       }
      }
      }
     else {                   #if y is input, then return the rectangular array
        for (i in 1: sizex ) {
            for (j in 1:sizey) {
                 result[j,i]<-sum((!is.na(x[,i]))&(!is.na(y[,j])))
                     }
                } }
   return (result) }
   
   
  #minor routines to manipulate matrices for factor analysis 
#output from eigen is two lists: $values and $vectors
#
eigen.loadings <- function (x) { #convert eigen vectors to loadings by unnormalizing them
    n <- length(x$values)
    x$values[ x$values<0]<- 0
    fix<-sqrt(x$values)
    result<- x$vectors * rep(fix, each = n)}
    
    
 #do a principal components analysis by doing eigen vector decomposition and then showing the biggest nfactors 
principal <-
function(r,nfactors=0) {
   n <- dim(r)[1]
    result<-list(values=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),residual=matrix(rep(0,n*n),ncol=n),fit=0)
    eigens <-eigen(r)    #call the eigen value decomposition routine
    result$values <-eigens$values
    result$loadings <-eigen.loadings(eigens)
    rownames(result$loadings) <- rownames(r)
    if (nfactors>0) result$loadings<-result$loadings[,1:nfactors]
    residual<- factor.residuals(r,result$loadings)
    r2<- sum(r*r)
    rstar2<- sum(residual*residual)
    result$residual<-residual
    result$fit<- 1-rstar2/r2
    return(result)
   }
   
   
   
   
 #function to find factor congruence ( 
 
 factor.congruence <- function (x,y) {
  nx<- dim(x)[2]
   ny<- dim(y)[2]
   cross<- t(y) %*% x   #inner product will have dim of ny * nx
   sumsx<- sqrt(1/diag(t(x)%*%x))   
   sumsy<- sqrt(1/diag(t(y)%*%y)) 

   result<- matrix(rep(0,nx*ny),ncol=nx)
    result<-  sumsy * (cross * rep(sumsx, each = ny))
   return(result)
   }
   
   
 #estimate a correlation matrix given a particular factor/component model 
 factor.model <- function(x) { 
    result<- x%*%t(x)
    return (result)}
  
  #find residuals
  factor.residuals <- function(r, f) {
   rstar<- r- factor.model(f)
   return(rstar)}
   
   factor.fit <-function (r,f) {
     r2 <-sum( r*r)
     rstar <- factor.residuals(r,f)
     rstar2 <- sum(rstar*rstar)
     fit<- 1- rstar2/r2
     return(fit) }
     
 
  #do a principal components analysis by doing eigen vector decomposition and then showing the biggest nfactors 
 principal <- function(r,nfactors=0) {
   n <- dim(r)[1]
    result<-list(values=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),residual=matrix(rep(0,n*n),ncol=n),fit=0)
    eigens <-eigen(r)    #call the eigen value decomposition routine
    result$values <-eigens$values
    result$loadings <-eigen.loadings(eigens)
    if (nfactors>0) result$loadings<-result$loadings[,1:nfactors]
    residual<- factor.residuals(r,result$loadings)
    r2<- sum(r*r)
    rstar2<- sum(residual*residual)
    result$residual<-residual
    result$fit<- 1-rstar2/r2
    return(result)
   }
   
   
   factor.rotate <-function(f,angle,col1,col2)  {
     nvar<- dim(f)[2]
     rot<- matrix(rep(0,nvar*nvar),ncol=nvar)
     rot[cbind(1:nvar, 1:nvar)] <- 1
     theta<- 2*pi*angle/360
     rot[col1,col1]<-cos(theta)
     rot[col2,col2]<-cos(theta)
     rot[col1,col2]<- -sin(theta)
     rot[col2,col1]<- sin(theta)
     result <- f %*% rot
     return(result) }
  
#phi coefficient from a 2 x 2 table  (function probably already exists but I can't find it)
phi <- function(t) {  #t is a 2 x 2 matrix
   
    r.sum<-rowSums(t)
    c.sum <-colSums(t)
    total<-sum(r.sum)
    r.sum<- r.sum/total
    c.sum <-c.sum/total
    v<- prod(r.sum,c.sum)
    phi <- (t[1,1]/total- c.sum[1]*r.sum[1])/sqrt(v)
  }
    
   #an alternative function to do the same, need to combine into one 
  phi1 <- function(x) {  # returns the phi coefficient   x  is a set of 4 numbers e.g. c(1,2,34)  
   xm <- matrix(x,ncol=2)
   total <- sum(xm)
   xp <- xm/total
   HR <- xp[1,1] + xp[1,2]
   SR <- xp[1,1]+ xp[2,1]
   VP <- xp[1,1]
   phi <- (VP - HR*SR) /sqrt(HR*(1-HR)*(SR)*(1-SR))
   return(phi)
   }
  
 #convert a  phi coefficient to the corresponding polychoric correlation 
 #needs John Fox's polychor program
 phi2poly <- function(ph,cp,cc) {
     #ph is the phi coefficient
     #cp is the selection ratio of the predictor
     #cc is the success rate of the criterion
     r.marg<-rep(0,2)
     c.marg<- rep(0,2)
     p<-array(rep(0,4),dim=c(2,2))
     r.marg[1]<- cp
     r.marg[2]<- 1 -cp 
     c.marg[1]<- cc
     c.marg[2]<- 1-cc
     
     p[1,1]<- r.marg[1]*c.marg[1]+ ph*sqrt(prod(r.marg,c.marg))
     p[2,2]<- r.marg[2]*c.marg[2]+ ph*sqrt(prod(r.marg,c.marg))
     p[1,2]<- r.marg[1]*c.marg[2]- ph*sqrt(prod(r.marg,c.marg))
     p[2,1]<- r.marg[2]*c.marg[1]- ph*sqrt(prod(r.marg,c.marg))
     
     result<-polychor(p ) 
     return(result)}
     
     
 psycho.demo <- function() {    
 #simulate correlation matrix with variable cut points -- psychometric demo
 #make up some correlations with different cut points
cuts <-c(-2,-1,0,1,2)
nsub<-1000    #how many subjects
r<-.6         #population correlation of observed with theta
theta <-rnorm(nsub)  #make up some random normal theta scores
err<- rnorm(nsub)    #random normal error scores

obser<- theta*(r) + err*sqrt(1-r*r)  #observed = T + E

#convert to 0/1  with different values of cuts
trunc<- matrix(rep(obser,length(cuts)),ncol=length(cuts))  
for (i in 1:length(cuts)) {
   trunc[obser>cuts[i],i]<- 1
   trunc[obser< cuts[i],i]<- 0}
   
d.mat<- data.frame(theta,obser,trunc)  #combine into a data frame
trunc.cor<- cor(d.mat)                 #find the correlations
freq <- mean(d.mat)                    #find the frequencies of scores

#now, convert the upper diagonal to polychorics using John Fox's polychor and my phi2poly

for (i in 4:length(d.mat)) {
   for (j in 3:i) {
       trunc.cor[j,i]<- phi2poly(trunc.cor[i,j],freq[i],freq[j]) 
       }}

} 


#quick demonstration of a random walk process 
#to show range of variability over trials and replications
#the "confidence lines" represent 2 sd of theoretical sum 

random.walk <- function(length=200,n=20,effect=0) {
	num <- seq(1:length)
	colors <- rainbow(n)                 #select colors
 	x <- cumsum(rnorm(length,effect))/num
 	plot(x,ylim = c(-2.5,2.5),typ="b",col=colors[1],main="Sample means as function of sample size",ylab="sample mean",xlab="sample size")
	for (i in 2:n) {
 		x <- cumsum(rnorm(length,effect))/num     #find the next line
 		points(x,col=colors[i],typ="b")         #draw it
 			}
  	curve(2/sqrt(x),1,length,101,add=TRUE)     #upper confidence region
 	curve(-2/sqrt(x),1,length,101,add=TRUE)    #lower confidence region
 	text(length/2,2,"effect size is " )
 	text(length/2, 1.8,eval(effect))
 }

#read from clipboard for Mac (thanks to Ken Knoblauch for this hint
#for PCs, use read.table(file("clipboard"))
read.clipboard<-function(header=TRUE,...) {
    MAC<-Sys.info()[1]=="Darwin"    #are we on a Mac using the Darwin system?
   if (!MAC ) {if (header) read.clipboard<-read.table(file("clipboard"),header=TRUE,...)
            else read.clipboard<-read.table(file("clipboard"),...) }
    else {
   if (header) read.clipboard<-  read.table(pipe("pbpaste"),header=TRUE,...)
   else read.clipboard<- read.table(pipe("pbpaste") ,...)}
   }
   
 read.clipboard.csv<-function(header=TRUE,sep=',',...) {  #same as read.clipboard(sep=',') 
    MAC<-Sys.info()[1]=="Darwin"    #are we on a Mac using the Darwin system?
   if (!MAC ) {if (header) read.clipboard<-read.table(file("clipboard"),header=TRUE,...)
            else read.clipboard<-read.table(file("clipboard"),...) }
    else {
   if (header) read.clipboard<-  read.table(pipe("pbpaste"),header=TRUE,...)
   else read.clipboard<- read.table(pipe("pbpaste") ,...)}
   }



#functions to help analyze data from correlation matrices
 #used as part of the Synthetic Aperture Personality Assessment  (SAPA) project
 #
 #extract cluster keys from a factor/principal components loadings matrix
 
 factor2cluster <- 
function(r,nfactors=0) {
   n <- dim(r)[1]
    result<-list(values=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),residual=matrix(rep(0,n*n),ncol=n),fit=0)
    eigens <-eigen(r)    #call the eigen value decomposition routine
    result$values <-eigens$values
    result$loadings <-eigen.loadings(eigens)
    rownames(result$loadings) <- rownames(r)
    if (nfactors>0) result$loadings<-result$loadings[,1:nfactors]
    residual<- factor.residuals(r,result$loadings)
    r2<- sum(r*r)
    rstar2<- sum(residual*residual)
    result$residual<-residual
    result$fit<- 1-rstar2/r2
    return(result)
   }

#find the correlation matrix of scales made up of items defined in a keys matrix (e.g., extracted by factor2cluster) 
#takes as input the keys matrix as well as a correlation matrix of all the items
cluster.cor <- 
    function(keys,r.mat,correct=TRUE) { #function to extract clusters according to the key vector
				#default is to correct for attenuation and show this above the diagonal
				#find the correlation matrix of scales made up of items defined in a keys matrix (e.g., extracted by factor2cluster) 
                #takes as input the keys matrix as well as a correlation matrix of all the items
 if(!is.matrix(keys)) keys <- as.matrix(keys)  #keys are sometimes a data frame - must be a matrix
 covar <- t(keys) %*% r.mat %*% keys    #matrix algebra is our friend
 var <- diag(covar)
 sd.inv <- 1/sqrt(var)
 ident.sd <- diag(sd.inv,ncol = length(sd.inv))
 cluster.correl <- ident.sd %*% covar  %*% ident.sd
 key.var <- diag(t(keys) %*% keys)
 key.alpha <- ((var-key.var)/var)*(key.var/(key.var-1))
 key.alpha[is.nan(key.alpha)] <- 1           #if only 1 variable to the cluster, then alpha is undefined
 key.alpha[!is.finite(key.alpha)] <- 1   
 colnames(cluster.correl) <- names(key.alpha)
 rownames(cluster.correl) <- names(key.alpha)
 if (correct) {cluster.corrected <- correct.cor(cluster.correl,t(key.alpha))}  #correct for attenuation
 
 return(list(cor=cluster.correl,sd=sqrt(var),corrected= cluster.corrected))
 }
 
 
 
 #a function to extract subsets of variables (a and b) from a correlation matrix m
     #and find the multiple correlation beta weights + R2 of the a set predicting the b set
     mat.regress <- function(m,a,b) {
        #first reorder the matrix to select the right variables
         nm <- dim(m)[1]
        t.mat <- matrix(0,ncol=nm,nrow=nm)
        ab <- c(a,b)
         numa <- length(a)
     	numb <- length(b)
        nab <- numa+numb
        for (i in 1:nab) {
     	t.mat[i,ab[i]] <- 1 }
     	
     	reorder <- t.mat %*% m %*% t(t.mat)
     	a.matrix <- reorder[1:numa,1:numa]
     	b.matrix <- reorder[1:numa,(numa+1):nab]
     	model.mat <- solve(a.matrix,b.matrix)       #solve the equation bY~aX
     	if (length(b) >1 ) { rownames(model.mat) <- rownames(m)[a]
     	 colnames(model.mat) <- colnames(m)[b]
     	 
     	r2 <- colSums(model.mat * b.matrix) }
     	 else { r2 <- sum(model.mat * b.matrix)
     	 names(model.mat) <- rownames(m)[a]
     	 names(r2) <- colnames(m)[b]}
     	mat.regress <- list(beta=model.mat,r2=r2)
     	return(mat.regress)
     	}
