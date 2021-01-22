library(S4)
library(mclust)#mclust has a function for calculating ARI
######################################################################
#Estimation K for Kmeans for Settings A1-A7 and B1-B5. For all 11 methods including S4
#Number of estimating K, ARI, RMSE of K and average of estimating K are shown
#Reproduce the Table 2 and Web Table 5 in the paper
######################################################################
num_cores=50
num_simulation=100
n.resample=100

temp.k<-matrix(0,nrow=11,ncol=num_simulation)
temp.ari<-matrix(0,nrow=11,ncol=num_simulation)

RMSE<-function(x,y){
  return(sqrt(sum((x-y)^2)/length(x)))
}
res.k<-matrix(0,nrow=12,ncol=11)
res.ari<-matrix(0,nrow=12,ncol=11)
res.avek<-matrix(0,nrow=12,ncol=11)
res.rmsek<-matrix(0,nrow=12,ncol=11)
numk<-c(3,4,4,4,4,2,2,4,4,4,2,2)#true number of clusters for each setting

setting.name<-c(paste("A",1:7,sep=""),paste("B",1:5,sep=""))
for(ii in setting.name){
  wcs = mclapply(1:num_simulation,function(i){
    set.seed(i)
    X<-Sim1(settings = ii)
    x1<-X$x
    label<-X$label
    
    #library(mclust)
    gapst_pca <- K.Clust(x1,method="Gap.PCA",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
    gapst_pca<-gapst_pca$model
    gapst_pca$Tab<-gapst_pca$Tab[-1,]
    gaprs_pca <- (2:10)[with(gapst_pca,maxSE(Tab[,"gap"],Tab[,"SE.sim"]))]
    res<-kmeans(x1,centers=gaprs_pca,nstart=50)
    ari_gaprs_pca<-adjustedRandIndex(label,res$cluster)#adjustedRandIndex is from package mclust
    
    
    gaprs_unif <- K.Clust(x1,method="Gap.Unif",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
    gaprs_unif<-gaprs_unif$model
    gaprs_unif$Tab<-gaprs_unif$Tab[-1,]
    gaprs_unif <- (2:10)[with(gaprs_unif,maxSE(Tab[,"gap"],Tab[,"SE.sim"]))]
    res<-kmeans(x1,centers=gaprs_unif,nstart=50)
    ari_gaprs_unif<-adjustedRandIndex(label,res$cluster)
    
    res.PS <- K.Clust(x1,method="PS",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
    res.PS<-res.PS$model
    ps<-(2:10)[max(which(res.PS$mean.pred[-1]==max(res.PS$mean.pred[-1])))]
    res<-kmeans(x1,centers=ps,nstart=50)
    ari_ps<-adjustedRandIndex(label,res$cluster)
    
    
    jm = K.Clust(x1,method="Jump",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
    jm<-jm$model
    jm<-(2:10)[which.max(jm$jumps[-1])]
    res<-kmeans(x1,centers=jm,nstart=50)
    ari_jm<-adjustedRandIndex(label,res$cluster)
    
    CH<-K.Clust(x1,method="CH",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
    res<-kmeans(x1,centers=CH,nstart=50)
    ari_CH<-adjustedRandIndex(label,res$cluster)
    
    rboot<-K.Clust(x1,method="FW",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
    rboot<-rboot$model
    rboot<-(2:10)[max(which(rboot$stabk[-1]==max(rboot$stabk[-1])))]
    res<-kmeans(x1,centers=rboot,nstart=50)
    ari_rboot<-adjustedRandIndex(label,res$cluster)
    
    
    res.LD<-K.Clust(x1,method="LD",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
    res<-kmeans(x1,centers=res.LD,nstart=50)
    ari_LD<-adjustedRandIndex(label,res$cluster)
    
    
    
    res.S4 <- K.Clust(x1,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
    res<-kmeans(x1,centers=res.S4$maxk,nstart=50)
    ari_S4<-adjustedRandIndex(label,res$cluster)
    
    
    KL<-K.Clust(x1,method="KL",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
    res<-kmeans(x1,centers=KL,nstart=50)
    ari_KL<-adjustedRandIndex(label,res$cluster)
    
    
    H<-K.Clust(x1,method="H",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
    res<-kmeans(x1,centers=H,nstart=50)
    ari_H<-adjustedRandIndex(label,res$cluster)
    
    
    
    Sil<-K.Clust(x1,method="silhouette",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
    res<-kmeans(x1,centers=Sil,nstart=50)
    ari_Sil<-adjustedRandIndex(label,res$cluster)
    
    res.ari<-c(ari_gaprs_pca,ari_gaprs_unif,ari_jm,
               ari_CH,ari_KL,ari_H,ari_Sil,ari_ps,ari_rboot,ari_LD,ari_S4)
    res.k<-c(gaprs_pca,gaprs_unif,jm,CH,KL,H,Sil,ps,rboot,res.LD,res.S4$maxk)
    
    return(list(res.ari=res.ari,res.k=res.k))
    
  },mc.cores = num_cores)
  for(i in 1:num_simulation){
    temp.k[,i]<-wcs[[i]]$res.k
    temp.ari[,i]<-wcs[[i]]$res.ari
  }
  res.k[which(ii==setting.name),]<-apply(temp.k,1,function(x){sum(x==numk[which(ii==setting.name)])})
  res.ari[which(ii==setting.name),]<-apply(temp.ari,1,mean)
  res.avek[which(ii==setting.name),]<-apply(temp.k,1,mean)
  res.rmsek[which(ii==setting.name),]<-apply(temp.k,1,function(x){RMSE(x,numk[which(ii==setting.name)])})
}



colnames(res.k)<-c("gap_pca","gap_unif","jump","CH","KL","H","Sil",
                   "PS","FW","LD","S4")
colnames(res.ari)<-c("gap_pca","gap_unif","jump","CH","KL","H","Sil",
                     "PS","FW","LD","S4")
colnames(res.avek)<-c("gap_pca","gap_unif","jump","CH","KL","H","Sil",
                      "PS","FW","LD","S4")
colnames(res.rmsek)<-c("gap_pca","gap_unif","jump","CH","KL","H","Sil",
                       "PS","FW","LD","S4")

write.csv(res.k,file="res.k.csv")
write.csv(res.ari,file="res.ari.csv")
write.csv(res.avek,file="res.avek.csv")
write.csv(res.rmsek,file="res.rmsek.csv")
