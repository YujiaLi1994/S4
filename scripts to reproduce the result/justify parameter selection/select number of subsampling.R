num_cores=50
num_simulation=100
num.B<-100



numk<-c(3,4,4,4,4,2,2,4,4,4,2,2)

#----------------------------------B=20
n.resample=20
wcs20<-c()
ARI_list20<-list()
setting.name<-c(paste("A",1:7,sep=""),paste("B",1:5,sep=""))
for(ii in setting.name){
  
  ARI<-matrix(NA,nrow=num_simulation,ncol=num.B)
  wcs = mclapply(1:num_simulation,function(i){
    set.seed(i)
    X<-Sim1(settings = ii)
    x1<-X$x
    label<-X$label
    ari<-c()
    for(ind.B in 1:num.B){
      set.seed(ind.B)
      res<-K.Clust(x1,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
      model<-kmeans(x1,centers=res$maxk,nstart=20)
      ari[ind.B]<-adjustedRandIndex(model$cluster,label)
    }
    
    return(ari)
  },mc.cores = num_cores)
  for(i in 1:num_simulation){
    ARI[i,]<-wcs[[i]]
  }
  ARI_list20[[which(ii==setting.name)]]<-ARI
  #wcs20[as.numeric(ii)]<-sum(unlist(wcs==numk[as.numeric(ii)]))
  wcs20[which(ii==setting.name)]<-mean(apply(ARI,1,var))
  
}

# wcs20_new<-c()
# for(i in 1:12){
#   ARI<-ARI_list20[[i]]
#   wcs20_new[i]<-mean(apply(round(ARI,5),1,function(x){1-sum((table(x)/length(x))^2)}))
# }

write.csv(wcs20,file="S4_Select_subsampling_B20.csv")

#----------------------------------B=50
n.resample=50
wcs50<-c()
ARI_list50<-list()
for(ii in setting.name){
  
  ARI<-matrix(NA,nrow=num_simulation,ncol=num.B)
  wcs = mclapply(1:num_simulation,function(i){
    set.seed(i)
    X<-Sim1(settings = ii)
    x1<-X$x
    label<-X$label
    ari<-c()
    for(ind.B in 1:num.B){
      set.seed(ind.B)
      res<-K.Clust(x1,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
      model<-kmeans(x1,centers=res$maxk,nstart=20)
      ari[ind.B]<-adjustedRandIndex(model$cluster,label)
    }
    
    return(ari)
  },mc.cores = num_cores)
  for(i in 1:num_simulation){
    ARI[i,]<-wcs[[i]]
  }
  ARI_list50[[which(ii==setting.name)]]<-ARI
  #wcs20[as.numeric(ii)]<-sum(unlist(wcs==numk[as.numeric(ii)]))
  wcs50[which(ii==setting.name)]<-mean(apply(ARI,1,var))
  
}

# wcs50_new<-c()
# for(i in 1:12){
#   ARI<-ARI_list50[[i]]
#   wcs50_new[i]<-mean(apply(round(ARI,5),1,function(x){1-sum((table(x)/length(x))^2)}))
# }

write.csv(wcs50,file="S4_Select_subsampling_B50.csv")

#----------------------------------B=100
n.resample=100
wcs100<-c()
ARI_list100<-list()
for(ii in setting.name){
  
  ARI<-matrix(NA,nrow=num_simulation,ncol=num.B)
  wcs = mclapply(1:num_simulation,function(i){
    set.seed(i)
    X<-Sim1(settings = ii)
    x1<-X$x
    label<-X$label
    ari<-c()
    for(ind.B in 1:num.B){
      set.seed(ind.B)
      res<-K.Clust(x1,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
      model<-kmeans(x1,centers=res$maxk,nstart=20)
      ari[ind.B]<-adjustedRandIndex(model$cluster,label)
    }
    
    return(ari)
  },mc.cores = num_cores)
  for(i in 1:num_simulation){
    ARI[i,]<-wcs[[i]]
  }
  ARI_list100[[which(ii==setting.name)]]<-ARI
  #wcs20[as.numeric(ii)]<-sum(unlist(wcs==numk[as.numeric(ii)]))
  wcs100[which(ii==setting.name)]<-mean(apply(ARI,1,var))
  
}

# wcs100_new<-c()
# for(i in 1:12){
#   ARI<-ARI_list100[[i]]
#   wcs100_new[i]<-mean(apply(round(ARI,5),1,function(x){1-sum((table(x)/length(x))^2)}))
# }

write.csv(wcs100,file="S4_Select_subsampling_B100.csv")

#----------------------------------B=200
n.resample=200
wcs200<-c()
ARI_list200<-list()
for(ii in setting.name){
  
  ARI<-matrix(NA,nrow=num_simulation,ncol=num.B)
  wcs = mclapply(1:num_simulation,function(i){
    set.seed(i)
    X<-Sim1(settings = ii)
    x1<-X$x
    label<-X$label
    ari<-c()
    for(ind.B in 1:num.B){
      set.seed(ind.B)
      res<-K.Clust(x1,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
      model<-kmeans(x1,centers=res$maxk,nstart=20)
      ari[ind.B]<-adjustedRandIndex(model$cluster,label)
    }
    
    return(ari)
  },mc.cores = num_cores)
  for(i in 1:num_simulation){
    ARI[i,]<-wcs[[i]]
  }
  ARI_list200[[which(ii==setting.name)]]<-ARI
  #wcs20[as.numeric(ii)]<-sum(unlist(wcs==numk[as.numeric(ii)]))
  wcs200[which(ii==setting.name)]<-mean(apply(ARI,1,var))
  
}

# wcs200_new<-rep(NA,12)
# for(i in 1:11){
#   ARI<-ARI_list200[[i]]
#   wcs200_new[i]<-mean(apply(round(ARI,5),1,function(x){1-sum((table(x)/length(x))^2)}))
# }
write.csv(wcs200,file="S4_Select_subsampling_B200.csv")

#----------------------------------B=500
n.resample=500
wcs500<-c()
ARI_list500<-list()
for(ii in setting.name){
  
  ARI<-matrix(NA,nrow=num_simulation,ncol=num.B)
  wcs = mclapply(1:num_simulation,function(i){
    set.seed(i)
    X<-Sim1(settings = ii)
    x1<-X$x
    label<-X$label
    ari<-c()
    for(ind.B in 1:num.B){
      set.seed(ind.B)
      res<-K.Clust(x1,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
      model<-kmeans(x1,centers=res$maxk,nstart=20)
      ari[ind.B]<-adjustedRandIndex(model$cluster,label)
    }
    
    return(ari)
  },mc.cores = num_cores)
  for(i in 1:num_simulation){
    ARI[i,]<-wcs[[i]]
  }
  ARI_list500[[which(ii==setting.name)]]<-ARI
  #wcs20[as.numeric(ii)]<-sum(unlist(wcs==numk[as.numeric(ii)]))
  wcs500[which(ii==setting.name)]<-mean(apply(ARI,1,var))
  
}
# wcs500_new<-c()
# for(i in 1:12){
#   ARI<-ARI_list500[[i]]
#   wcs500_new[i]<-mean(apply(round(ARI,5),1,function(x){1-sum((table(x)/length(x))^2)}))
# }

write.csv(wcs500,file="S4_Select_subsampling_B500.csv")

save(ARI_list20,file="ARI_list20.Rdata")
save(ARI_list50,file="ARI_list50.Rdata")
save(ARI_list100,file="ARI_list100.Rdata")
save(ARI_list200,file="ARI_list200.Rdata")
save(ARI_list500,file="ARI_list500.Rdata")

#wcs_new<-cbind(wcs20_new,wcs50_new,wcs100_new,wcs200_new)
wcs<-cbind(wcs20,wcs50,wcs100,wcs200,wcs500)
