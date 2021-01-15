##' Estimate number of clusters and sparsity parameter simultaneously by extended gap statistic.
##' @title Estimate number of clusters and sparsity parameter simultaneously by extended gap statistic.
##' @param x data matrix, where rows represent samples and columns represents features
##' @param lambda_list a list whose length is equal to k_vector. Every element is a vector of lambda
##' as the search region for this K.
##' @param k_vector The search pool for number of clusters
##' @param n.perms Number of permutated data generated
##' @param num.cores Number of cpus used in parallel computing. It is suggest to run large num_cores
##' in servers instead of laptops.
##' @return A list of two components:
##' \itemize{
##' \item{res_k: }{Optimal number of clusters}
##' \item{res_w: }{Optimal lambda}
##'
##' }
##' @references Estimating the number of clusters in a data set via the gap statistic. Journal of the Royal Statistical Society: Series B (Statistical Methodology),63(2):411-423.
##'
##' Witten, D. M. and Tibshirani, R. (2010). A framework for feature selection in clustering. Journal of the American Statistical Association, 105(490):713-726.
##'
##' @export
##' @examples
##' \dontrun{
##'#generate simulation II data
##'data<-Sim2(h=200,q = 50,u=0.8)#For demo purpose, we use 200 features in total for fast result.
##'#using Efficient Algorithm for Choosing Grids of lambda for each K
##'k_vector<-2:7#search K from 2 to 7
##'wbounds_list<-list(1)
##'for(l in 1:length(k_vector)){#for each K, using the algorithm to get 20 lambda.
##'  wbounds_list[[l]] = region.lambda(lam1=1.5,iteration=20,data,k_vector[l])
##'}

##'#This part of code is to delete the lambda which selecting all the genes.
##'#since our S4 method calculate specificity, the lambda selecting all the genes must be removed.
##'for(l in 1:length(k_vector)){
##'  temp<-KMeansSparseCluster(data,K=k_vector[l],wbounds=wbounds_list[[l]],nstart=100)
##'  num<-rep(0,length(temp))
##'  for(i in 1:length(num)){
##'    num[i]<-sum(temp[[i]]$ws>0)
##'  }
##'  if(sum(num==ncol(data))>0){
##'    wbounds_list[[l]]<-wbounds_list[[l]][1:(min(which(num==ncol(data)))-3)]
##'  }
##'}

##'#run extended Gap statistic method
##'res.Gap<-KL.Gap(x=data,lambda_list = wbounds_list,k_vector = k_vector,n.perm = 50,num.cores = 1)
##' }





KL.Gap<-function(x,k_vector=c(2,3,4),lambda_list,n.perms=25,num.cores=1){
  res_gap<-list(1)
  res_gap_sd<-list(1)
  res = mclapply(1:length(k_vector),function(l){
    a = KmeansSparseCluster.permute1(x, K=k_vector[l], nperms = n.perms, wbounds = lambda_list[[l]],
                                     silent = FALSE, nvals = 20, centers=NULL)
    return(list(gap=a$gaps,sd=a$sdgaps))
  },mc.cores = num.cores)


  for(i in 1:length(res)){
    res_gap[[i]]<-res[[i]]$gap
    res_gap_sd[[i]]<-res[[i]]$sd
  }


  maxgap_k<-rep(-1,length(k_vector))
  for(i in 1:length(k_vector)){
    maxgap_k[i]<-max(res_gap[[i]])
  }
  gap_k<-k_vector[max(which(max(maxgap_k)==maxgap_k))]
  #a<-KmeansSparseCluster.permute1(x, K=gap_k, nperms = n.perms, wbounds = lambda_list[[which(k_vector==gap_k)]],
  #                                silent = FALSE, nvals = 20, centers=NULL)
  #maxgap_minus_sd<-max(a$gaps)-a$sdgaps[which.max(a$gaps)]
  a<-res_gap[[which(k_vector==gap_k)]]
  maxgap_minus_sd<-max(a)-res_gap_sd[[which(k_vector==gap_k)]][which.max(a)]
  res_w<-lambda_list[[which(k_vector==gap_k)]][which.max(a>maxgap_minus_sd)]
  return(list(res_k=gap_k,res_w=res_w))
}
