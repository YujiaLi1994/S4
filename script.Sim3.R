library(KLClust)
num_simulation=10
num_cores=10
n.resample = 25 #resampling times
prop = 0.7 #resample n proportion
trim=0.05
n.perms<-25
n.div<-25




#------------------------------------------------------
#1 ari
#2 jaccard
#3 choose K
#4 number of features selected

gapwcs1<-matrix(rep(-1,2*num_simulation),ncol=num_simulation)
gapwcs2<-matrix(rep(-1,2*num_simulation),ncol=num_simulation)
gapwcs3<-matrix(rep(-1,2*num_simulation),ncol=num_simulation)
gapwcs4<-matrix(rep(-1,2*num_simulation),ncol=num_simulation)
s4wcs1<-matrix(rep(-1,2*num_simulation),ncol=num_simulation)
s4wcs2<-matrix(rep(-1,2*num_simulation),ncol=num_simulation)
s4wcs3<-matrix(rep(-1,2*num_simulation),ncol=num_simulation)
s4wcs4<-matrix(rep(-1,2*num_simulation),ncol=num_simulation)
pswcs1<-matrix(rep(-1,2*num_simulation),ncol=num_simulation)
pswcs2<-matrix(rep(-1,2*num_simulation),ncol=num_simulation)
pswcs3<-matrix(rep(-1,2*num_simulation),ncol=num_simulation)
pswcs4<-matrix(rep(-1,2*num_simulation),ncol=num_simulation)

# rwcs1<-matrix(rep(-1,4*4),ncol=4)
# rwcs2<-matrix(rep(-1,4*4),ncol=4)
# rwcs3<-matrix(rep(-1,4*4),ncol=4)
# rwcs4<-matrix(rep(-1,4*4),ncol=4)

POIN1<-40
POIN2<-30
POIN3<-20
sigma1<-0.2
sigma2<-1
POIM<-20
Ua<-4
Ub<-10
j1<-1
j2<-3
Nnoise<-100



#k_vector<-c(2,3,4,5,6,7)
k_vector<-c(3)
#for (j in 1:5){

for (k in 1:2){
  for(ii in 1:num_simulation){
    q = c(10) #number of DE features out of 1000
    cov_list<-c(0.1,0.3,0.5)
    #########cov=0.1
    if(j2==1){
      ulow = c(0.8,1)
      uupper=c(1.0,1.5)
    }
    if(j2==2){
      ulow = c(1.5,2)
      uupper=c(2,2.5)
    }
    ##########cov=0.3
    
    ###############cov=0.5
    if(j2==3){
      ulow = c(1.5,2)
      uupper=c(2,2.5)
    }
    #h<-1000#DE evidence
    K<-3
    
    set.seed(ii)
    
    data<-Sim3(POIN1=POIN1,POIN2=POIN2,POIN3=POIN3,Nm=POIM,
               M=q[j1],Ua=Ua,Ub=Ub,Ulow=ulow[k],Uupper=uupper[k],sigma1=sigma1,sigma2=sigma2,Nnoise=Nnoise,cov=cov_list[j2])
    x<-data$data
    wbounds_list<-list(1)
    for(l in 1:length(k_vector)){
      wbounds_list[[l]] = region.lambda(lam1=3,iteration=10,x=x,K=k_vector[l])
    }
    #wbounds = region.lamda(b1=6,b2=5,x=x,h=h,K=3)
    ##test whether it is dense enough
    for(l in 1:length(k_vector)){
      temp<-KMeansSparseCluster(x,K=k_vector[l],wbounds=wbounds_list[[l]],nstart=100)
      num<-rep(0,length(temp))
      for(i in 1:length(num)){
        num[i]<-sum(temp[[i]]$ws>0)
      }
      if(sum(num==ncol(x))>0){
        wbounds_list[[l]]<-wbounds_list[[l]][1:(min(which(num==ncol(x)))-3)]
      }
    }
    
    ###############gap 
    ##number of feature selected
    res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms = n.perms,num.cores=num_cores)
    
    result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                     silent = FALSE, maxiter=6, centers=NULL)
    
    #compute ARI
    true.label = data$true_label #underlying truth
    clus_label_gap = result.gap[[1]]$Cs # class label
    ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
    
    TrueFeature<-data$true_feature
    SelectedFeature_gap <- as.integer(result.gap[[1]]$ws>0)
    
    #compute jaccard index
    jaccard_gap<-Jaccard.index(TrueFeature,SelectedFeature_gap)
    gap_k<-res.gap$res_k
    num_gap<-sum(result.gap[[1]]$ws>0)####number of feature selected
    
    ################S4
    # Estimate wbounds using my method, K fixed at 3
    res.s4 = KL.S4(x, lambda_list=wbounds_list,trim=0.05,k_vector=k_vector,n.resample=n.resample,num.cores=num_cores)
    result.s4 = KMeansSparseCluster(x, K=res.s4$optimal_k, wbounds = res.s4$optimal_lambda, nstart = 100, ######
                                     silent = FALSE, maxiter=6, centers=NULL)
    
    #compute ARI
    true.label = data$true_label #underlying truth
    clus_label_s4 = result.s4[[1]]$Cs # class label
    ari_s4<-adjustedRand(true.label,clus_label_s4)[[2]]
    
    TrueFeature<-data$true_feature
    #TrueFeature[1:q[j]] <- 1
    SelectedFeature_s4 <- as.integer(result.s4[[1]]$ws>0)
    
    #compute jaccard index
    jaccard_s4<-Jaccard.index(TrueFeature,SelectedFeature_s4)
    s4_k<-res.s4$optimal_k
    num_s4<-sum(result.s4[[1]]$ws>0)####number of feature selected
    ########for prediction strength
    res.ps = KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num_cores=num_cores)
    
    result.ps = KMeansSparseCluster(x, K=res.ps$res_k, wbounds = res.ps$res_w, nstart = 100, ######
                                    silent = FALSE, maxiter=6, centers=NULL)
    
    #compute ARI
    true.label = data$true_label #underlying truth
    clus_label_ps = result.ps[[1]]$Cs # class label
    ari_ps<-adjustedRand(true.label,clus_label_ps)[[2]]
    
    TrueFeature<-data$true_feature
    
    SelectedFeature_ps <- as.integer(result.ps[[1]]$ws>0)
    
    #compute jaccard index
    jaccard_ps<-Jaccard.index(TrueFeature,SelectedFeature_ps)
    ps_k<-res.ps$res_k
    num_ps<-sum(result.ps[[1]]$ws>0)####number of feature selected
    
    s4wcs1[k,ii]<-ari_s4
    s4wcs2[k,ii]<-jaccard_s4
    s4wcs3[k,ii]<-s4_k
    s4wcs4[k,ii]<-num_s4
    gapwcs1[k,ii]<-ari_gap
    gapwcs2[k,ii]<-jaccard_gap
    gapwcs3[k,ii]<-gap_k
    gapwcs4[k,ii]<-num_gap
    pswcs1[k,ii]<-ari_gap
    pswcs2[k,ii]<-jaccard_gap
    pswcs3[k,ii]<-gap_k
    pswcs4[k,ii]<-num_gap

    
    
    
  }
}