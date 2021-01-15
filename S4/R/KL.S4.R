##' Estimate number of clusters and sparsity parameter simultaneously by S4.
##' @title Estimate number of clusters and sparsity parameter simultaneously by S4.
##' @param x data matrix, where rows represent samples and columns represents features
##' @param lambda_list a list whose length is equal to k_vector. Every element is a vector of lambda
##' as the search region for this K.
##' @param k_vector The search pool for number of clusters
##' @param trim Percentage of trim by S4 method.
##' @param n.resample Number of subsampling, Usually needed to be larger than 20
##' @param num.cores Number of cpus used in parallel computing. It is suggest to run large num_cores
##' in servers instead of laptops.
##' @return A list of two components:
##' \itemize{
##' \item{clus_score: }{Cluster stability score for every K and lambda.}
##' \item{final_score: }{The average of cluster stability score and feature stability
##' score for every K and lambda}
##' \item{optimal_k: }{Optimal number of clusters}
##' \item{optimal_lambda: }{Optimal lambda}
##' }
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

##'#run S4 method
##'res.S4<-KL.S4(x=data,lambda_list = wbounds_list,k_vector = k_vector,trim =0.05,n.resample = 50,num.cores = 1)
##' }




KL.S4 = function(x, lambda_list,trim=0.05,k_vector=c(2,3,4),n.resample=25,num.cores=1){
  result_clus<-list(1)
  result_final<-list(1)
  for(ind_k in 1:length(k_vector)){
    k<-k_vector[ind_k]
    lambda_vector<-lambda_list[[ind_k]]
    wcs = mclapply(1:length(lambda_vector),function(ind_lambda,k,x){
      lambda<-lambda_vector[ind_lambda]
      error = T
      ind<-0
      while(error){
        result= tryCatch(score_lambda(lambda,x,k,trim=0.05,n.resample=n.resample), error = function(x) NA)
        if(!is.na(result[1])){
          error = F
        }
        ind<-ind+1
        if(ind>3){break}
      }


      #result<-score_lambda(lambda,x,k,trim=0.05,n.resample=25)
      return(result)
    },k = k,x=x,mc.cores = num.cores)
    result_clus[[ind_k]]<-NA
    result_final[[ind_k]]<-NA
    for(ind_lambda in 1:length(lambda_vector)){
      result_clus[[ind_k]][ind_lambda]<-wcs[[ind_lambda]]$clus_score
      result_final[[ind_k]][ind_lambda]<-wcs[[ind_lambda]]$final_score
    }
  }
  result_k<-rep(-1,length(k_vector))
  for(i in 1:length(k_vector)){
    result_k[i]<-max(result_clus[[i]])
  }
  res_k<-k_vector[max(which(result_k==max(result_k)))]
  index<-which.max(result_final[[which(res_k==k_vector)]])
  res_lambda<-lambda_list[[which(res_k==k_vector)]][index]
  return(list(clus_score=result_clus,final_score=result_final,optimal_k=res_k,optimal_lambda=res_lambda))
}
