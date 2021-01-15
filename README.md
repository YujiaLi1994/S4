# S4

This package includes the following parts. First, Three methods (extended gap statistic, extedned prediction strength and S4) for estimating `K` and `lambda` simultaneously in sparse K-means are included in `KL.Gap`, `KL.PS` and `KL.S4` respectively. Second, 11 methods for estimating K in K-means are included in `K.Clust` function. Third, 9 high-dimensional datasets are also included with name starting with `ds`. Several simulation functions are included in `Sim1`, `Sim2` and `Sim3`. The extensive simulation settings and real applications in this paper can also serve as a standard benchmark for evaluation purpose when new methods for simultaneous estimation of K and lambda are developed. Finally, an efficient algorithm for generate lambda grid for each K is included in the `region.lambda`.

The R scripts contains all the code to reproduce the result in the following paper:
Paper: Li, Yujia, et al. "Simultaneous Estimation of Number of Clusters and Feature Sparsity Parameter in High-dimensional clustering analysis". Biometrics, accepted (2021). A previous arxiv version: arXiv:1909.01930

## Installing
S4 package files are in `S4/` folder, You can install by copying and paste the following code into R

```
devtools::install_github("YujiaLi1994/S4/S4")
```
Alternatively, download the `tar.gz` zipped file, unzip it and install ogClust package using `devtools::install()`, which can make sure to install all the dependences. Make sure R packages mvtnorm, caret, parallel, irr, plyr,cluster, fpc, MCMCpack, MASS, sparcl have all been properly imported.


## Estimating Number of Cluster K for K-means
```r
#simulate the data of simulation I by Sim1 function.
data<-Sim1(settings = "A1")#A1-A7 are well separated data and B1-B7 are not-well separated data
#using S4
#trim.S4 is the percentage of trimmed scattered points.
#Cutoff is the threshold to detect the K=1
#n.resample is number of subsampling
res.S4<-K.Clust(data$x,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=50)
#users can also implement other methods, see the document and reference for details.
#Take prediction strength for an example
res.PS<-K.Clust(data$x,method="PS",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=50)
``` 
## Estimating Number of Clusters K and sparsity parameter lambda simultaneously for sparse K-means
*Simulate a dataset using Sim2 function.

```r
data<-Sim2(h=200,q = 50,u=0.8)#For demo purpose, we use 200 features in total for fast result. 

k_vector<-2:7#K is selected from 2 to 7
```
* Use Efficient Algorithm for Choosing Grids of lambda for each K
```r
wbounds_list<-list(1)
for(l in 1:length(k_vector)){#for each K, using the algorithm to get 20 lambda.
  wbounds_list[[l]] = region.lambda(lam1=1.5,iteration=20,data,k_vector[l])
}
```
* Delete lambda which will select all the feature, since our S4 method calculate specificity, the lambda selecting all the genes must be removed. Including very large lambda (that select all feature) has the possibility to cause an error.
```r
for(l in 1:length(k_vector)){
  temp<-KMeansSparseCluster(data,K=k_vector[l],wbounds=wbounds_list[[l]],nstart=100)
  num<-rep(0,length(temp))
  for(i in 1:length(num)){#get the corresponding number of features for each K and each lambda
    num[i]<-sum(temp[[i]]$ws>0)
  }
  if(sum(num==ncol(data))>0){#For each K, if a certain lambda selects all features, delete it. 
  #For a large simulation study, two closest lambda next to it can also be removed to be conservative.
    wbounds_list[[l]]<-wbounds_list[[l]][1:(min(which(num==ncol(data)))-3)]
  }
}
```
* implementing S4
```r
#run S4 method
res.S4<-KL.S4(x=data,lambda_list = wbounds_list,k_vector = k_vector,trim =0.05,n.resample = 50,num.cores = 1)
```

* Implement extended Gap and extended Prediction strength
```r
#run extended Gap statistic method
res.Gap<-KL.Gap(x=data,lambda_list = wbounds_list,k_vector = k_vector,n.perm = 50,num.cores = 1)
#run extended Prediction strength method
res.PS<-KL.PS(x=data,lambda_list = wbounds_list,k_vector = k_vector,cv = 2,M=20,num.cores = 1,cutoff = 0.8)
```
##Access nine real datasets used in S4 paper
All the datasets are imported and accessible after library(S4).
```r
#Data sets in package ‘S4’:
# ds.GSE13159                  Microarray data of Leukemia after preprocessing(GSE13159)
# ds.GSE17855                  Microarray data of leukemia after preprocessing(GSE17855)
# ds.GSE47474                  RNA Sequencing data of rat brain after
# preprocessing(GSE47474)
# ds.GSE6891                   Microarray data of leukemia after preprocessing(GSE6891)
# ds.ISOLET                    letter-name dataset after preprocessing
# ds.Leaf                      Plant species leaves dataset after preprocessing
# ds.Pancancer                 RNA Sequencing data after preprocessing(Pancancer)
# ds.SNP                       SNP data after preprocessing
# ds.TissueType                Microarray data after preprocessing(Mammalian tissue
#                                                                  types dataset:)
```
