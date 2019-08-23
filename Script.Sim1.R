library(KLClust)
##############test low dimensional function
num_cores=10
num_simulation=100
n.resample=25
wcs1 = mclapply(1:num_simulation,function(i){
  x1<-Sim1(settings = "1")
  
  gaprs_pca <- K.Clust(x1,method="Gap.PCA",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  gaprs_unif <- K.Clust(x1,method="Gap.Unif",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  ps <- K.Clust(x1,method="PS",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  jm = K.Clust(x1,method="Jump",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  CH<-K.Clust(x1,method="CH",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  rboot<-K.Clust(x1,method="FW",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  re_sens<-K.Clust(x1,method="LD",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  #taylor_null <- taylor.clustering_trim(x1,K.try = 2:10,n.resample=n.resample,specificity=TRUE,n.null=n.null,subnull=TRUE,trim=trim)
  res.S4 <- K.Clust(x1,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  KL<-K.Clust(x1,method="KL",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  H<-K.Clust(x1,method="H",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  Sil<-K.Clust(x1,method="silhouette",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  return(list(gap_pca=gaprs_pca,gap_unif=gaprs_unif,jump=jm,CH=CH,KL=KL,H=H,Sil=Sil,
              prediction_strength=ps,rboot=rboot,re_sens=re_sens,S4=res.S4))
  
},mc.cores = num_cores)


wcs2 = mclapply(1:num_simulation,function(i){
  x1<-Sim1(settings = "2")
  
  gaprs_pca <- K.Clust(x1,method="Gap.PCA",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  gaprs_unif <- K.Clust(x1,method="Gap.Unif",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  ps <- K.Clust(x1,method="PS",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  jm = K.Clust(x1,method="Jump",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  CH<-K.Clust(x1,method="CH",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  rboot<-K.Clust(x1,method="FW",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  re_sens<-K.Clust(x1,method="LD",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  #taylor_null <- taylor.clustering_trim(x1,K.try = 2:10,n.resample=n.resample,specificity=TRUE,n.null=n.null,subnull=TRUE,trim=trim)
  res.S4 <- K.Clust(x1,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  KL<-K.Clust(x1,method="KL",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  H<-K.Clust(x1,method="H",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  Sil<-K.Clust(x1,method="silhouette",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  return(list(gap_pca=gaprs_pca,gap_unif=gaprs_unif,jump=jm,CH=CH,KL=KL,H=H,Sil=Sil,
              prediction_strength=ps,rboot=rboot,re_sens=re_sens,S4=res.S4))
  
},mc.cores = num_cores)

wcs3 = mclapply(1:num_simulation,function(i){
  x1<-Sim1(settings = "3")
  
  gaprs_pca <- K.Clust(x1,method="Gap.PCA",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  gaprs_unif <- K.Clust(x1,method="Gap.Unif",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  ps <- K.Clust(x1,method="PS",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  jm = K.Clust(x1,method="Jump",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  CH<-K.Clust(x1,method="CH",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  rboot<-K.Clust(x1,method="FW",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  re_sens<-K.Clust(x1,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  #taylor_null <- taylor.clustering_trim(x1,K.try = 2:10,n.resample=n.resample,specificity=TRUE,n.null=n.null,subnull=TRUE,trim=trim)
  res.S4 <- K.Clust(x1,method="KL",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  KL<-K.Clust(x1,method="KL",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  H<-K.Clust(x1,method="H",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  Sil<-K.Clust(x1,method="silhouette",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  return(list(gap_pca=gaprs_pca,gap_unif=gaprs_unif,jump=jm,CH=CH,KL=KL,H=H,Sil=Sil,
              prediction_strength=ps,rboot=rboot,re_sens=re_sens,S4=res.S4))
  
},mc.cores = num_cores)


wcs4 = mclapply(1:num_simulation,function(i){
  x1<-Sim1(settings = "4")
  
  gaprs_pca <- K.Clust(x1,method="Gap.PCA",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  gaprs_unif <- K.Clust(x1,method="Gap.Unif",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  ps <- K.Clust(x1,method="PS",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  jm = K.Clust(x1,method="Jump",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  CH<-K.Clust(x1,method="CH",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  rboot<-K.Clust(x1,method="FW",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  re_sens<-K.Clust(x1,method="LD",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  #taylor_null <- taylor.clustering_trim(x1,K.try = 2:10,n.resample=n.resample,specificity=TRUE,n.null=n.null,subnull=TRUE,trim=trim)
  res.S4 <- K.Clust(x1,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  KL<-K.Clust(x1,method="KL",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  H<-K.Clust(x1,method="H",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  Sil<-K.Clust(x1,method="silhouette",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  return(list(gap_pca=gaprs_pca,gap_unif=gaprs_unif,jump=jm,CH=CH,KL=KL,H=H,Sil=Sil,
              prediction_strength=ps,rboot=rboot,re_sens=re_sens,S4=res.S4))
  
},mc.cores = num_cores)

wcs5 = mclapply(1:num_simulation,function(i){
  x1<-Sim1(settings = "5")
  
  gaprs_pca <- K.Clust(x1,method="Gap.PCA",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  gaprs_unif <- K.Clust(x1,method="Gap.Unif",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  ps <- K.Clust(x1,method="PS",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  jm = K.Clust(x1,method="Jump",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  CH<-K.Clust(x1,method="CH",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  rboot<-K.Clust(x1,method="FW",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  re_sens<-K.Clust(x1,method="LD",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  #taylor_null <- taylor.clustering_trim(x1,K.try = 2:10,n.resample=n.resample,specificity=TRUE,n.null=n.null,subnull=TRUE,trim=trim)
  res.S4 <- K.Clust(x1,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  KL<-K.Clust(x1,method="KL",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  H<-K.Clust(x1,method="H",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  Sil<-K.Clust(x1,method="silhouette",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  return(list(gap_pca=gaprs_pca,gap_unif=gaprs_unif,jump=jm,CH=CH,KL=KL,H=H,Sil=Sil,
              prediction_strength=ps,rboot=rboot,re_sens=re_sens,S4=res.S4))
  
},mc.cores = num_cores)


wcs6 = mclapply(1:num_simulation,function(i){
  x1<-Sim1(settings = "6")
  
  gaprs_pca <- K.Clust(x1,method="Gap.PCA",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  gaprs_unif <- K.Clust(x1,method="Gap.Unif",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  ps <- K.Clust(x1,method="PS",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  jm = K.Clust(x1,method="Jump",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  CH<-K.Clust(x1,method="CH",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  rboot<-K.Clust(x1,method="FW",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  re_sens<-K.Clust(x1,method="LD",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  #taylor_null <- taylor.clustering_trim(x1,K.try = 2:10,n.resample=n.resample,specificity=TRUE,n.null=n.null,subnull=TRUE,trim=trim)
  res.S4 <- K.Clust(x1,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  KL<-K.Clust(x1,method="KL",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  H<-K.Clust(x1,method="H",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  Sil<-K.Clust(x1,method="silhouette",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  return(list(gap_pca=gaprs_pca,gap_unif=gaprs_unif,jump=jm,CH=CH,KL=KL,H=H,Sil=Sil,
              prediction_strength=ps,rboot=rboot,re_sens=re_sens,S4=res.S4))
  
},mc.cores = num_cores)


wcs7 = mclapply(1:num_simulation,function(i){
  x1<-Sim1(settings = "7")
  
  gaprs_pca <- K.Clust(x1,method="Gap.PCA",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  gaprs_unif <- K.Clust(x1,method="Gap.Unif",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  ps <- K.Clust(x1,method="PS",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  jm = K.Clust(x1,method="Jump",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  CH<-K.Clust(x1,method="CH",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  rboot<-K.Clust(x1,method="FW",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  re_sens<-K.Clust(x1,method="LD",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  #taylor_null <- taylor.clustering_trim(x1,K.try = 2:10,n.resample=n.resample,specificity=TRUE,n.null=n.null,subnull=TRUE,trim=trim)
  res.S4 <- K.Clust(x1,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  KL<-K.Clust(x1,method="KL",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  H<-K.Clust(x1,method="H",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  Sil<-K.Clust(x1,method="silhouette",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  return(list(gap_pca=gaprs_pca,gap_unif=gaprs_unif,jump=jm,CH=CH,KL=KL,H=H,Sil=Sil,
              prediction_strength=ps,rboot=rboot,re_sens=re_sens,S4=res.S4))
  
},mc.cores = num_cores)


wcs8 = mclapply(1:num_simulation,function(i){
  x1<-Sim1(settings = "8")
  
  gaprs_pca <- K.Clust(x1,method="Gap.PCA",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  gaprs_unif <- K.Clust(x1,method="Gap.Unif",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  ps <- K.Clust(x1,method="PS",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  jm = K.Clust(x1,method="Jump",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  CH<-K.Clust(x1,method="CH",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  rboot<-K.Clust(x1,method="FW",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  re_sens<-K.Clust(x1,method="LD",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  #taylor_null <- taylor.clustering_trim(x1,K.try = 2:10,n.resample=n.resample,specificity=TRUE,n.null=n.null,subnull=TRUE,trim=trim)
  res.S4 <- K.Clust(x1,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  KL<-K.Clust(x1,method="KL",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  H<-K.Clust(x1,method="H",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  Sil<-K.Clust(x1,method="silhouette",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  return(list(gap_pca=gaprs_pca,gap_unif=gaprs_unif,jump=jm,CH=CH,KL=KL,H=H,Sil=Sil,
              prediction_strength=ps,rboot=rboot,re_sens=re_sens,S4=res.S4))
  
},mc.cores = num_cores)


wcs9 = mclapply(1:num_simulation,function(i){
  x1<-Sim1(settings = "9")
  
  gaprs_pca <- K.Clust(x1,method="Gap.PCA",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  gaprs_unif <- K.Clust(x1,method="Gap.Unif",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  ps <- K.Clust(x1,method="PS",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  jm = K.Clust(x1,method="Jump",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  CH<-K.Clust(x1,method="CH",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  rboot<-K.Clust(x1,method="FW",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  re_sens<-K.Clust(x1,method="LD",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  #taylor_null <- taylor.clustering_trim(x1,K.try = 2:10,n.resample=n.resample,specificity=TRUE,n.null=n.null,subnull=TRUE,trim=trim)
  res.S4 <- K.Clust(x1,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  KL<-K.Clust(x1,method="KL",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  H<-K.Clust(x1,method="H",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  Sil<-K.Clust(x1,method="silhouette",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  return(list(gap_pca=gaprs_pca,gap_unif=gaprs_unif,jump=jm,CH=CH,KL=KL,H=H,Sil=Sil,
              prediction_strength=ps,rboot=rboot,re_sens=re_sens,S4=res.S4))
  
},mc.cores = num_cores)


wcs10 = mclapply(1:num_simulation,function(i){
  x1<-Sim1(settings = "10")
  
  gaprs_pca <- K.Clust(x1,method="Gap.PCA",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  gaprs_unif <- K.Clust(x1,method="Gap.Unif",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  ps <- K.Clust(x1,method="PS",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  jm = K.Clust(x1,method="Jump",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  CH<-K.Clust(x1,method="CH",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  rboot<-K.Clust(x1,method="FW",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  re_sens<-K.Clust(x1,method="LD",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  #taylor_null <- taylor.clustering_trim(x1,K.try = 2:10,n.resample=n.resample,specificity=TRUE,n.null=n.null,subnull=TRUE,trim=trim)
  res.S4 <- K.Clust(x1,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  KL<-K.Clust(x1,method="KL",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  H<-K.Clust(x1,method="H",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  Sil<-K.Clust(x1,method="silhouette",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  return(list(gap_pca=gaprs_pca,gap_unif=gaprs_unif,jump=jm,CH=CH,KL=KL,H=H,Sil=Sil,
              prediction_strength=ps,rboot=rboot,re_sens=re_sens,S4=res.S4))
  
},mc.cores = num_cores)





nc = c("1","2","3","4","5","6","7","8","9","10")

awc1 = matrix(unlist(wc1),ncol = Number_of_Simulation)
awc2 = matrix(unlist(wc2),ncol = Number_of_Simulation)
awc3 = matrix(unlist(wc3),ncol = Number_of_Simulation)
awc4 = matrix(unlist(wc4),ncol = Number_of_Simulation)
awc5 = matrix(unlist(wc5),ncol = Number_of_Simulation)
awc6 = matrix(unlist(wc6),ncol = Number_of_Simulation)
awc7 = matrix(unlist(wc7),ncol = Number_of_Simulation)
awc8 = matrix(unlist(wc8),ncol = Number_of_Simulation)
awc9 = matrix(unlist(wc9),ncol = Number_of_Simulation)
awc10 = matrix(unlist(wc10),ncol = Number_of_Simulation)
# awc11 = matrix(unlist(wc11),ncol = Number_of_Simulation)
# awc12 = matrix(unlist(wc12),ncol = Number_of_Simulation)
# 
# 
# 
# #wc18 = matrix(unlist(wc8),ncol = Number_of_Simulation)
# 
# 
# #wc = rbind(wc11,wc12,wc13,wc14,wc15,wc16,wc17,wc18)
wc = rbind(awc1,awc2,awc3,awc4,awc5,awc6,awc7,awc8,awc9,awc10)
#wc = wc17
#wc = rbind(wc12,wc13,wc14,wc15,wc16,wc17)
wc = t(apply(wc,1,table))

for (i in 1:length(wc)){
  wc[[i]][nc[!nc%in%names(wc[[i]])]] <- 0
  wc[[i]] <- wc[[i]][nc]
}

wc = matrix(unlist(wc),ncol = 10,byrow = T)

rownames(wc) = rep(c("Gap_pca","gap_unif","jump","CH","KL","H","Sil","Pred str","FW","LD","S4"),10)
#rownames(wc) = rep(c("Gap","Pred str","jump","per","unif","pca"),8)
wc
#write.csv(wc,file="100simulationlowdimension12315.csv")