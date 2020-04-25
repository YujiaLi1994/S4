# Project Title

This package includes the following parts. 
First, Three methods (gap statistic, prediction strength and S4) for estimating K and lambda simultaneously in sparse K-means are included in KL.Gap, KL.PS and KL.S4 respectively. 
Second, 11 methods for estimating K in K-means are included in K.Clust function. 
Third, nine high-dimensional datasets are also included with name starting with ds. 
Finally, an efficient algorithm for generate Lambda grid for each K is included in the region.lambda 
and several simulation functions are included in Sim1, Sim2 and Sim3.

The R script contains all the code to reproduce the result in the following paper:
Paper: Li, Yujia, et al. "Simultaneous Estimation of Number of Clusters and Feature Sparsity in Clustering High-Dimensional Data." arXiv preprint arXiv:1909.01930 (2019).

### Installing

For installation, we recommend to unzip the tar.gz file first and then use devtools::install() to install the package, which can make sure to install all the depends. Make sure R package mvtnorm, caret, clues, parallel, irr, plyr,cluster, fpc, MCMCpack, MASS, sparcl, ggplot2 have all been properly imported.

## Below is the code example:

## Estimating number of Cluster in K-means
  library(S4)
  
  x1<-Sim1(settings = "1")  
  
  gaprs_pca <- K.Clust(x1,method="Gap.PCA",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  gaprs_unif <- K.Clust(x1,method="Gap.Unif",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  ps <- K.Clust(x1,method="PS",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  jm = K.Clust(x1,method="Jump",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  CH<-K.Clust(x1,method="CH",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  rboot<-K.Clust(x1,method="FW",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  re_sens<-K.Clust(x1,method="LD",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  res.S4 <- K.Clust(x1,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  KL<-K.Clust(x1,method="KL",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  H<-K.Clust(x1,method="H",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  Sil<-K.Clust(x1,method="silhouette",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
## Estimating Cluster number and sparsity parameter simultaneously for sparse K-means
library(S4)

n.resample=50

n.div=50

n.perms=50

x<-ds.GSE17855$data

k_vector<-c(2,3,4,5,6,7)

## get lambda tuning list
wbounds_list<-mclapply(1:length(k_vector),function(i){
  K<-k_vector[i]
  error = T
  ind<-0
  while(error){
    result= tryCatch(region.lambda(lam1=3,iteration=30,x,K), error = function(x) NA)
    if(!is.na(result[1])){
      error = F
    }
    ind<-ind+1
    if(ind>3){break}
  }
  return(result)
  
},mc.cores = 20)

for(l in 1:length(k_vector)){
  temp<-KMeansSparseCluster(x,K=k_vector[l],wbounds=wbounds_list[[l]],nstart=100)
  num<-rep(0,length(temp))
  for(i in 1:length(num)){
    num[i]<-sum(temp[[i]]$ws>0)
  }
  if(sum(num==ncol(x))>0){
    wbounds_list[[l]]<-wbounds_list[[l]][1:(min(which(num==ncol(x)))-2)]
  }
}
## implementing S4, extended Gap and extended Prediction strength
res.s4 = KL.S4(x, lambda_list=wbounds_list,trim=0.05,k_vector=k_vector,n.resample=n.resample,num.cores=20)

res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num.cores=10)

res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num.cores=10)
