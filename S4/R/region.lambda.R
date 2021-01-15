##' Efficient Algorithm for Choosing Grids of lambda
##' @title Efficient Algorithm for Choosing Grids of lambda
##' @param lam1 The lower bound for lambda
##' @param iteration Number of iterative search 
##' @param x data matrix, where rows represent samples and columns represents features 
##' @param K Number of clusters


##' @return A vector of sorted lambda value

##' @export
##' @examples
##' \dontrun{
##'data(data_tanbin_mat)
##'data<-RatBrain$data
##'region.lambda(lam1=3,iteration=20,x=data,K=3)
##' }






region.lambda<-function(lam1=3,iteration=30,x,K){
  lam1<-lam1
  h<-ncol(x)
  lam2<-sqrt(h)
  lam<-(lam1+lam2)/2
  iter<-2
  lam.vector<-c(lam1,lam2)
  n.vector<-rep(-1,length(lam.vector))
  for(i in 1:length(n.vector)){
    temp1<-KMeansSparseCluster1(x, K=K, nstart=100,wbounds = lam.vector[i])##############
    n.vector[i]<-sum(temp1[[1]]$ws!=0)
  }
  n.vector_old<-n.vector
  lam<-(lam.vector[1]+lam.vector[2])/2
  while(iter<=iteration){
    print(iter)
    lam.vector<-sort(c(lam.vector,lam))
    n.vector<-rep(-1,length(lam.vector))
    n.vector[-which(lam.vector==lam)]<-n.vector_old
    temp1<-KMeansSparseCluster1(x, K=K, nstart=100,wbounds = lam.vector[which(lam.vector==lam)])##############
    n.vector[which(lam.vector==lam)]<-sum(temp1[[1]]$ws!=0)
    d<-rep(-1,(length(lam.vector)-1))
    for(i in 1:length(d)){
      if(n.vector[i]==n.vector[i+1]){
        d[i]<-0
      }else{
        d[i]<-log2(n.vector[i+1]/n.vector[i])
      }
      
    }
    index<-which.max(d)
    lam<-(lam.vector[index]+lam.vector[index+1])/2
    iter<-iter+1
    n.vector_old<-n.vector
  }
  lam.vector<-sort(c(lam.vector,lam))
  
  return(lam.vector)
}
