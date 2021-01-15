##' Estimate number of clusters for K-Means
##' @title Estimate number of clusters for K-Means
##' @param data A matrix where rows represent samples and columns represent features
##' @param method Method to use, including all the methods compared in S4 paper. Please refer to the reference of S4 paper
##' \itemize{
##' \item{"S4": }{Simultaneous Estimation of Number of Clusters and Feature Sparsity parameter in High-Dimensional Data Clustering Analysis. Biometrics accepted}
##' \item{"Gap.PCA", "Gap.Unif": }{Tibshirani, R., Walther, G., and Hastie, T. (2001).}
##' \item{"PS": }{Tibshirani, R. and Walther, G. (2005). }
##' \item{"Jump": }{Sugar, C. A. and James, G. M. (2003)}
##' \item{"CH": }{Calinski, T. and Harabasz, J. (1974)}
##' \item{"FW": }{Fang, Y. and Wang, J. (2012)}
##' \item{"LD": }{Levine, E. and Domany, E. (2001).}
##' \item{"KL": }{Krzanowski, W. J. and Lai, Y. (1988).}
##' \item{"H": }{Hartigan, J. A. (1975). }
##' \item{"silhouette": }{Rousseeuw, P. J. (1987).}
##' }
##' @param Kmin integer. Minimum number of clusters
##' @param Kmax integer. Maximum number of clusters.
##' @param trim.S4 The precentage of trim for S4 method
##' @param cutoff The threshold for detecting null cluster structure for resampling-based methods. If the maximum score for all K is smaller than cutoff, return one cluster.
##' @param n.resample Number of resamples in resampling-based methods. Also for Gap statistic, n.resample means the number of reference data to generate.

##' @return For S4 method, "maxk" is the optimal number of clusters and "score" is the stability score for K from kmin to kmax. For Gap, PS, Jump and FW, "maxk" is the estimated number of clusters and model is the original output (These methods are implemented internally using other packages (package "fpc" and package "cluster", readers can also read the reference for details). For other methods, only optimal number of clusters is returned.
##' @references S4: Simultaneous Estimation of Number of Clusters and Feature Sparsity parameter in High-Dimensional Data Clustering Analysis. Biometrics accepted.
##'
##' "Gap.PCA" and "Gap.Unif": Tibshirani, R., Walther, G., and Hastie, T. (2001). Estimating the number of clusters in a data set via the gap statistic. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 63(2):411-423.
##'
##' "PS":  Tibshirani, R. and Walther, G. (2005). Cluster validation by prediction strength. Journal of Computational and Graphical Statistics, 14(3):511-528.
##'
##' "Jump":  Sugar, C. A. and James, G. M. (2003). Finding the number of clusters in a dataset: An informationtheoretic approach. Journal of the American Statistical Association, 98(463):750-763.
##'
##' "CH":  Calinski, T. and Harabasz, J. (1974). A dendrite method for cluster analysis. Communications in Statistics-theory and Methods, 3(1):1-27.
##'
##' "FW":  Fang, Y. and Wang, J. (2012). Selection of the number of clusters via the bootstrap method. Computational Statistics & Data Analysis, 56(3):468-477.
##'
##' "LD":  Levine, E. and Domany, E. (2001). Resampling method for unsupervised estimation of cluster validity. Neural computation, 13(11):2573{2593.
##'
##' "KL":  Krzanowski, W. J. and Lai, Y. (1988). A criterion for determining the number of groups in a dataset using sum-of-squares clustering. Biometrics, pages 23-34.
##'
##' "H":  Hartigan, J. A. (1975). Clustering algorithms, new york: John willey and sons. Inc. Pages113129.}
##'
##' "silhouette":  Rousseeuw, P. J. (1987). Silhouettes: a graphical aid to the interpretation and validation of cluster analysis. Journal of computational and applied mathematics, 20:53-65.
##' @export
##' @examples
##' \dontrun{
##' data<-Sim1(settings = "A1")
##' #using S4
##' res.S4<-K.Clust(data$x,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=50)
##' #Other methods Take prediction strength for an example.
##' res.PS<-K.Clust(data$x,method="PS",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=50)
##' }




K.Clust<-function(data,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=50){
  #method_list:(Gap.PCA, Gap.Unif, Jump, CH, KL, H, silhouette, PS, FW, LD, S4)
  #NBclust give more than 20 indice, For detail of index, please look at NBclust package, here we only consider 4 index method.
  if(method=="S4"){
    result.S4 <- S4(data,K.try = Kmin:Kmax,n.resample=n.resample,trim=trim.S4,cut.off=cutoff)
    return(result.S4)
  }
  if(method=="Gap.PCA"){
    gapst_pca <- clusGap(data, FUN = kmeans, nstart = 100, K.max = 10, B =n.resample,spaceH0 ="scaledPCA")
    gaprs_pca <- with(gapst_pca,maxSE(Tab[,"gap"],Tab[,"SE.sim"]))
    res<-list(maxk=gaprs_pca,model=gapst_pca)
    return(res)
  }
  if(method=="Gap.Unif"){
    gapst_Unif <- clusGap(data, FUN = kmeans, nstart = 100, K.max = Kmax, B = n.resample,spaceH0 ="original")
    gaprs_Unif <- with(gapst_Unif,maxSE(Tab[,"gap"],Tab[,"SE.sim"]))
    res<-list(maxk=gaprs_Unif,model=gapst_Unif)
    return(res)
  }
  if(method=="PS"){
    ps <- prediction.strength(data,Gmin=Kmin,Gmax=Kmax,M=n.resample,cutoff=cutoff)
    res<-list(maxk=ps$optimalk,model=ps)
    #ps <- ps$optimalk
    return(res)
  }
  if(method=="Jump"){
    jm = jump(data,K=Kmax)
    res<-list(maxk=jm$maxjump,model=jm)
    return(res)
  }
  if(method=="CH"){
    CH<-cluster_CH(data,K.try=Kmin:Kmax)
    return(CH)
  }
  if(method=="FW"){
    rboot<-nselectboot(data,B=n.resample,clustermethod=kmeansCBI,krange=Kmin:Kmax)
    res<-list(maxk=rboot$kopt,model=rboot)
    return(res)
  }
  if(method=="LD"){
    re_sens<-resample_sensitivity(data, K.try = Kmin:Kmax,n.resample=n.resample,cut.off = cutoff)
    return(re_sens)
  }
  if(method=="KL"){
    KL<-NBClust1(data=data,min.nc = Kmin, max.nc = Kmax,method="kmeans",index="kl")$Best.nc[1]
    return(KL)
  }
  if(method=="H"){
    H<-hartigan1(data,min.nc =Kmin,max.nc = Kmax)
    return(H)
  }
  if(method=="silhouette"){
    Sil<-NBClust1(data=data,min.nc = Kmin, max.nc = Kmax,method="kmeans",index="silhouette")$Best.nc[1]
    return(Sil)
  }


}
