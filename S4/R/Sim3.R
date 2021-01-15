##' Simulate high-dimensional data with correlated feature structure of three clusters (Simulation III).
##'
##' @title Simulate high-dimensional data with correlated feature structure of three clusters (Simulation III).
##' @param POIN1 Number of samples in cluster 1 follows possion distribution. POIN1 is possion mean.
##' @param POIN2 Number of samples in cluster 2 follows possion distribution. POIN2 is possion mean.
##' @param POIN3 Number of samples in cluster 3 follows possion distribution. POIN3 is possion mean.
##' @param M Number of gene modules of informative genes.
##' @param Nm Number of genes within each gene module follows possion distribution. Nm is possion mean.
##' @param Ua,Ub Template gene expression of each cluster for each informative gene is simulated from Unif(Ua,Ub)
##' @param Ulow,Uupper Effect size(The constraint for template gene expression).
##' The template gene expression of each two clusters for every gene will be between Ulow and Uupper, where Ulow is lower bound.
##' @param sigma1 Biological variation to template gene expression of informative genes
##' @param sigma2 Biological variation to template gene expression of noise genes
##' @param Nnoise Number of noise genes.
##' @param cov Determine the correlation of genes within each module. The higher cov, the larger correlation.
##' @return A list of three components:
##' \itemize{
##' \item{data: }{A data matrix where row represents samples and column represents features}
##' \item{true_label: }{class label}
##' \item{true_feature: }{a vector equal to number of columns of data, where 1 represents informative features and 0 otherwise.}
##' }
##' @references Manuscript: Simultaneous Estimation of Number of Clusters and Feature Sparsity in Clustering High-Dimensional Data Using Resampling.
##' @export
##' @examples
##' \dontrun{
##' #####this is the default hyperparameter in the simulation3 of manuscript.
##' POIN1<-40
##' POIN2<-30
##' POIN3<-20
##' sigma1<-0.2
##' sigma2<-1
##' POIM<-20
##' Ua<-4
##' Ub<-10
##' j1<-1
##' j2<-1
##' Nnoise<-600
##'
##' #################below are the tuning parameter
##' cov=0.1
##' Ulow=1.0
##' Uupper=1.5
##' M=10
##' Sim3(POIN1=POIN1,POIN2=POIN2,POIN3=POIN3,Nm=POIM,
##' M=M,Ua=Ua,Ub=Ub,Ulow=Ulow,Uupper=Uupper,sigma1=sigma1,sigma2=sigma2,Nnoise=Nnoise,cov=cov)
##' }








Sim3<-function(POIN1=POIN1,POIN2=POIN2,POIN3=POIN3,Nm=POIM,M,Ua,Ub,Ulow,Uupper,sigma1,sigma2,Nnoise=600,cov=0.5){
  ##sigma1 is the biological variation of predictive genes
  ##sigma2 is the bilogical variation of the noise genes
  K<-3
  N1<-rpois(1,POIN1)
  N2<-rpois(1,POIN2)
  N3<-rpois(1,POIN3)
  N<-N1+N2+N3

  indicator_sample<-rep(1:3,c(N1,N2,N3))
  num_genes_in_module<-rep(0,M)
  for(i in 1:length(num_genes_in_module)){
    num_genes_in_module[i]<-rpois(1,Nm)
    while(num_genes_in_module[i]==0){
      num_genes_in_module[i]<-rpois(1,Nm)
    }
  }
  #num_genes_in_module<-rpois(M,Nm)
  indicator_module<-rep(1:M,num_genes_in_module)
  template_expression<-matrix(rep(0,M*K),nrow=K)
  for(i in 1:M){
    template_expression[,i]<-runif(K,Ua,Ub)
    while(((max(template_expression[,i])-min(template_expression[,i]))<Ulow)|
          ((max(template_expression[,i])-min(template_expression[,i]))>Uupper)){
      template_expression[,i]<-runif(K,Ua,Ub)
    }

  }
  #apply(template_expression,2,function(x){return(max(x)-min(x))})
  #####add biological variation
  gene_expression<-matrix(rep(0,N*M),nrow=N)
  for(i in 1:M){
    for(j in 1:N){
      gene_expression[j,i]<-rnorm(1,template_expression[indicator_sample[j],i],sigma1)
    }
  }
  #apply(gene_expression[which(indicator_sample==1),],2,mean)
  data<-matrix(rep(0,N*length(indicator_module)),nrow=N)
  for(i in 1:N){
    for(m in 1:M){
      v1<-riwish(60,((1-cov)*diag(num_genes_in_module[m])+
                       cov*matrix(rep(1,num_genes_in_module[m]^2),ncol=num_genes_in_module[m])))
      v<-standardize_cov(v1)
      data[i,which(indicator_module==m)]<-rmvnorm(1,mean=rep(gene_expression[i,m],num_genes_in_module[m]),sigma=v)
    }
  }
  #######add housekeeping genes
  gene_expression_template_housekeeping<-runif(Nnoise,Ua,Ub)
  for(i in 1:Nnoise){
    data<-cbind(data,rnorm(nrow(data),mean=gene_expression_template_housekeeping[i],sigma2))
  }

  true_label<-c(rep(1,N1),rep(2,N2),rep(3,N3))
  true_feature<-c(rep(1,sum(num_genes_in_module)),rep(0,Nnoise))
  return(list(data=data,true_label=true_label,true_feature=true_feature))

}
