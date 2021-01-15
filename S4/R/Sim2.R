##' Simulate high-dimensional data with independent feature structure of three clusters. (Simulation II)
##' @title Simulate high-dimensional data with independent feature structure of three clusters.
##' @param h Total number of features
##' @param q Number of informative features
##' @param u Effect size
##' @return A matrix with 99 rows(samples) and h column(features),
##' Sample 1-33 is cluster 1, 34-66 is cluster 2 and 67-99 is cluster 3.
##' The first q features are informative fearures whereas others are random noise
##' @references Manuscript: Simultaneous Estimation of Number of Clusters and Feature Sparsity in Clustering High-Dimensional Data Using Resampling
##' @export
##' @examples
##' \dontrun{
##' h=1000
##' q=200
##' u=0.8
##' Sim2(h=h,q=200,u=0.8)
##' }






Sim2<-function(h,q,u){
  #q = c(50,200) #number of DE features out of 1000
  #u = c(1,0.8,0.6,0.4)
  #h<-1000#DE evidence

  x = rbind(matrix(rnorm(33*h),ncol = h),
            matrix(rnorm(33*h),ncol = h),
            matrix(rnorm(33*h),ncol = h))
  x[1:33,1:q] <- x[1:33,1:q]-u
  x[34:66,1:q] <- x[34:66,1:q]
  x[67:99,1:q] <- x[67:99,1:q]+u
  return(x)
}
