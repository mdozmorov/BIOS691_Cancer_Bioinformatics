### Complete Linkage: When merging clusters containing >1 item, use the maximum distance between items 
### as the distance between clusters, then find the clusters having the minimum distance.

### Set up the data
weight.kg<-c(15,49,13,45,85,66,12,10)
height.cm<-c(95,156,95,160,178,176,90,78)
data<-matrix(cbind(weight.kg,height.cm),nrow=8)
dimnames(data)[[2]]<-c("weight.kg","height.cm")

## Iteration 1, calculate distance between all singletons
dist<-dist(data,diag=TRUE)
round(dist,3) 
## Find the minimum distance 
min(as.dist(dist)[as.dist(dist)!=0])  # 2
### Merge 1 and 3 
### Iteration 2, update distances to the cluster that contains (1,3)
old.dist<-as.matrix(dist)
new.dist<-as.matrix(dist)
new.dist[1,3]<-new.dist[3,1]<-0
new.dist[1,2]<-new.dist[2,1]<-max(old.dist[1,2],old.dist[3,2])
new.dist[1,4]<-new.dist[4,1]<-max(old.dist[1,4],old.dist[3,4])
new.dist[1,5]<-new.dist[5,1]<-max(old.dist[1,5],old.dist[3,5])
new.dist[1,6]<-new.dist[6,1]<-max(old.dist[1,6],old.dist[3,6])
new.dist[1,7]<-new.dist[7,1]<-max(old.dist[1,7],old.dist[3,7])
new.dist[1,8]<-new.dist[8,1]<-max(old.dist[1,8],old.dist[3,8])
new.dist[3,2]<-new.dist[2,3]<-max(old.dist[1,2],old.dist[3,2])
new.dist[3,4]<-new.dist[4,3]<-max(old.dist[1,4],old.dist[3,4])
new.dist[3,5]<-new.dist[5,3]<-max(old.dist[1,5],old.dist[3,5])
new.dist[3,6]<-new.dist[6,3]<-max(old.dist[1,6],old.dist[3,6])
new.dist[3,7]<-new.dist[7,3]<-max(old.dist[1,7],old.dist[3,7])
new.dist[3,8]<-new.dist[8,3]<-max(old.dist[1,8],old.dist[3,8])
## Check that new distance matrix is symmetric
all.equal(new.dist,t(new.dist))
round(as.dist(new.dist),3)
min(as.dist(new.dist)[as.dist(new.dist)!=0])  # 5.656854

## Merge 2 and 4
## Iteration 3; Update distances to the cluster that contains 2 & 4
new.dist2<-new.dist
new.dist2[2,4]<-new.dist2[4,2]<-0
new.dist2[1,2]<-new.dist2[2,1]<-new.dist2[4,1]<-new.dist2[1,4]<-max(old.dist[1,2],old.dist[3,2],old.dist[1,4],old.dist[3,4])
new.dist2[3,2]<-new.dist2[2,3]<-new.dist2[4,3]<-new.dist2[3,4]<-max(old.dist[1,2],old.dist[3,2],old.dist[1,4],old.dist[3,4])
new.dist2[2,5]<-new.dist2[5,2]<-new.dist2[4,5]<-new.dist2[5,4]<-max(old.dist[2,5],old.dist[4,5])
new.dist2[2,6]<-new.dist2[6,2]<-new.dist2[4,6]<-new.dist2[6,4]<-max(old.dist[2,6],old.dist[4,6])
new.dist2[2,7]<-new.dist2[7,2]<-new.dist2[4,7]<-new.dist2[7,4]<-max(old.dist[2,7],old.dist[4,7])
new.dist2[2,8]<-new.dist2[8,2]<-new.dist2[4,8]<-new.dist2[8,4]<-max(old.dist[2,8],old.dist[4,8])
## Check that new distance matrix is symmetric
all.equal(new.dist2,t(new.dist2))
round(as.dist(new.dist2),3)
## Find the minimum distance for identifying what clusters are merged next
min(as.dist(new.dist2)[as.dist(new.dist2)!=0])  #5.830952

# Merge 7 with cluster 1,3
# Iteration 4, Update distances to cluster that contains 1,3,7
new.dist3<-new.dist2
new.dist3[1,7]<-new.dist3[7,1]<-new.dist3[3,7]<-new.dist3[7,3]<-0
new.dist3[1,2]<-new.dist3[2,1]<-new.dist3[4,1]<-new.dist3[1,4]<-max(old.dist[1,2],old.dist[3,2],old.dist[1,4],old.dist[3,4],old.dist[7,2],old.dist[7,4])
new.dist3[3,2]<-new.dist3[2,3]<-new.dist3[4,3]<-new.dist3[3,4]<-max(old.dist[1,2],old.dist[3,2],old.dist[1,4],old.dist[3,4],old.dist[7,2],old.dist[7,4])
new.dist3[7,2]<-new.dist3[2,7]<-new.dist3[4,7]<-new.dist3[7,4]<-max(old.dist[1,2],old.dist[3,2],old.dist[1,4],old.dist[3,4],old.dist[7,2],old.dist[7,4])
new.dist3[1,5]<-new.dist3[5,1]<-new.dist3[3,5]<-new.dist3[5,3]<-new.dist3[7,5]<-new.dist3[5,7]<-max(old.dist[1,5],old.dist[3,5],old.dist[7,5])
new.dist3[1,6]<-new.dist3[6,1]<-new.dist3[3,6]<-new.dist3[6,3]<-new.dist3[7,6]<-new.dist3[6,7]<-max(old.dist[1,6],old.dist[3,6],old.dist[7,6])
new.dist3[1,8]<-new.dist3[8,1]<-new.dist3[3,8]<-new.dist3[8,3]<-new.dist3[7,8]<-new.dist3[8,7]<-max(old.dist[1,8],old.dist[3,8],old.dist[7,8])
## Check that new distance matrix is symmetric
all.equal(new.dist3,t(new.dist3))
round(as.dist(new.dist3),3)
## Find the minimum distance for identifying what clusters are merged next
min(as.dist(new.dist3)[as.dist(new.dist3)!=0]) #17.72005
## Merge 8 with 1, 3, 7

## Now, rather than go through another iteration by hand, apply the hclust function and check results
hclust.complete<-hclust(dist(data),method="complete")
hclust.complete$height
# 2.000000   5.656854   5.830952  17.720045  19.104973  43.863424 125.000000
## Notice the minimum distances correspond to the complete linkage method, which is to use the max 
## distance between clusters as the measure of distance.
plot(hclust.complete)
## And the order of merging agrees with what we got by hand.