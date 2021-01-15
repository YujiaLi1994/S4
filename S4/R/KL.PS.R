##' Estimate number of clusters and sparsity parameter simultaneously by extended prediction strength.
##' @title Estimate number of clusters and sparsity parameter simultaneously by extended prediction strength.
##' @param x data matrix, where rows represent samples and columns represents features
##' @param lambda_list a list whose length is equal to k_vector. Every element is a vector of lambda
##' as the search region for this K.
##' @param k_vector The search pool for number of clusters
##' @param cv Number of cross validation to do in each resampling. cv=2 is suggested by
##' Tibshirani, R. and Walther, G. (2005)
##' @param M Number of resampling. Every resampling is to divide the data into train and half.
##' @param num.cores Number of cpus used in parallel computing. It is suggest to run large num_cores
##' in servers instead of laptops.
##' @return A list of two components:
##' \itemize{
##' \item{res_k: }{Optimal number of clusters}
##' \item{res_w: }{Optimal lambda}
##'
##' }
##' @references Tibshirani, R. and Walther, G. (2005). Cluster validation by prediction strength. Journal of Computational and Graphical Statistics, 14(3):511-528.
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

##'#run extended Prediction strength method
##'res.PS<-KL.PS(x=data,lambda_list = wbounds_list,k_vector = k_vector,cv = 2,M=20,num.cores = 1,cutoff=0.8)

##' }






KL.PS<-function(x,lambda_list,cv=2,k_vector=c(2,3,4),M=20,cutoff,num.cores=1){
  #n = nrow(x)
  result_k<-list()
  result_lambda<-list()
  for(i1 in 1:length(k_vector)){
    print(i1)
    K<-k_vector[i1]
    temp_k<-c()
    temp_lambda<-c()

    for(i2 in 1:length(wbounds_list[[i1]])){
      wbounds<-wbounds_list[[i1]][i2]
      res<-lapply(1:M,function(i1,cv,K,wbounds,x){
        error = T
        ind<-0
        while(error){
          result= tryCatch(score_PS(x,wbounds,K,cv), error = function(x) NA)
          result
          if(!is.na(result[1])){
            error = F
          }
          ind<-ind+1
          if(ind>3){break}
        }
        return(result)

      },cv=cv,K=K,wbounds=wbounds,x=x)
      res.final = sapply(res, function(x) x$final, simplify = "array")
      res.clus = sapply(res, function(x) x$clus, simplify = "array")
      temp_k[i2]<-mean(res.clus)
      temp_lambda[i2]<-mean(res.final)
      #temp<-score_PS(wbounds=wbounds_list[[i1]][i2],x=x,K=k,cv=2)

    }
    result_k[[i1]]<-temp_k
    result_lambda[[i1]]<-temp_lambda
  }

  score_k<-unlist(lapply(result_k,max))
  if(max(score_k<cutoff)){
    est_k<-k_vector[max(which(score_k==max(score_k)))]
  }else{
    est_k<-k_vector[max(which(score_k>cutoff))]
  }
  lambda<-lambda_list[[est_k-1]][which.max(result_lambda[[est_k-1]])]
  res<-list(res_k=est_k,res_w=lambda)
  return(res)
}
