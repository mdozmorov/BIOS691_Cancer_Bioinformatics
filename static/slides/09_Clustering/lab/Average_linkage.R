dist.avg<-as.matrix(dist(data, diag = TRUE, upper = TRUE))
dist.avg
new.dist<-matrix(rep(0,8*8),ncol=8)
for (i in c(2,4,5,6,7,8)) {
   for (j in c(1,3)) {
	new.dist[i,j]<-(dist.avg[i,1]+dist.avg[i,3])/2
	new.dist[j,i]<-(dist.avg[i,1]+dist.avg[i,3])/2
	}
	for (j in c(2,4,5,6,7,8) ) {
	new.dist[i,j]<-dist.avg[i,j]
	new.dist[j,i]<-dist.avg[j,i]
	}
}
new.dist
## new minimum is 5.464986  1,3,7

new.dist<-matrix(rep(0,8*8),ncol=8)
for (i in c(2,4,5,6,8)) {
   for (j in c(1,3,7)) {
	new.dist[i,j]<-(dist.avg[i,1]+dist.avg[i,3]+dist.avg[i,7])/2
	new.dist[j,i]<-(dist.avg[i,1]+dist.avg[i,3]+dist.avg[i,7])/2
	}
	for (j in c(2,4,5,6,8) ) {
	new.dist[i,j]<-dist.avg[i,j]
	new.dist[j,i]<-dist.avg[j,i]
	}
}
new.dist
## new minimum is 5.656854  2,4

