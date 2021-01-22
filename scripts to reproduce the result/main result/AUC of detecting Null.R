



num_cores=1
num_simulation=2
n.resample=20

library(MASS)
setting.name<-c(paste("A",1:7,sep=""),paste("B",1:5,sep=""))
final_res_S4<-matrix(0,12,num_simulation)
for(ii in setting.name){
  wcs = mclapply(1:num_simulation,function(i){
    set.seed(i)
    X<-Sim1(settings = ii)
    x1<-X$x
    res<-K.Clust(x1,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
    return(max(res$score))
  },mc.cores = num_cores)
  final_res_S4[which(ii==setting.name),]<-unlist(wcs)
}
#save(final_res_S4,file="S4_AUC_extra.Rdata")


Sim_Null<-function(dimension,option){
  if(option=="Unif"){
    x<-matrix(runif(200*dimension),ncol = dimension)
  }
  if(option=="Normal"){
    x<-matrix(rnorm(200*dimension),ncol = dimension)
  }
  return(x)
}
final_res_S4_Null<-matrix(0,12,num_simulation)
Sim_K_1<-function(dimension,option){
  wcs = mclapply(1:num_simulation,function(i){
    x<-Sim_Null(dimension=dimension,option=option)
    res<-K.Clust(x,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
    return(max(res$score))
  },mc.cores = num_cores)
  return(wcs)
}
final_res_S4_Null[1,]<-unlist(Sim_K_1(dimension=2,option="Unif"))
final_res_S4_Null[2,]<-unlist(Sim_K_1(dimension=3,option="Unif"))
final_res_S4_Null[3,]<-unlist(Sim_K_1(dimension=4,option="Unif"))
final_res_S4_Null[4,]<-unlist(Sim_K_1(dimension=5,option="Unif"))
final_res_S4_Null[5,]<-unlist(Sim_K_1(dimension=10,option="Unif"))
final_res_S4_Null[6,]<-unlist(Sim_K_1(dimension=20,option="Unif"))
final_res_S4_Null[7,]<-unlist(Sim_K_1(dimension=2,option="Normal"))
final_res_S4_Null[8,]<-unlist(Sim_K_1(dimension=3,option="Normal"))
final_res_S4_Null[9,]<-unlist(Sim_K_1(dimension=4,option="Normal"))
final_res_S4_Null[10,]<-unlist(Sim_K_1(dimension=5,option="Normal"))
final_res_S4_Null[11,]<-unlist(Sim_K_1(dimension=10,option="Normal"))
final_res_S4_Null[12,]<-unlist(Sim_K_1(dimension=20,option="Normal"))

#-----------------------------------------------calculate AUC
library(pROC)
AUC_S4<-matrix(0,12,12)
rownames(AUC_S4)<-c("unif_2","unif_3","unif_4","unif_5","unif_10","unif_20",
                    "normal_2","normal_3","normal_4","normal_5","normal_10","normal_20")
colnames(AUC_S4)<-setting.name
#pdf("S4.pdf")
par(mfrow=c(12,12))
for(i in 1:12){
  for(j in 1:12){
    S4_1<-as.numeric(final_res_S4[j,])
    S4_null<-as.numeric(final_res_S4_Null[i,])
    AUC_S4[i,j]<-auc(factor(c(rep("1",num_simulation),rep("0",num_simulation))),c(S4_1,S4_null),direction="<")
    #plot.roc(factor(c(rep("1",100),rep("0",100))),c(S4_1,S4_null),direction="<")
  }
  
}
#dev.off()
#--------------------Get AUC matrix for different dimension
#--------------------Web Table 10 shows the AUC using same dimensional nulls,
#--------------------user can extract the result from this table by the dimensions
write.csv(AUC_S4,file="S4_AUC.csv")

#----------------------------------Gap Unif
library(cluster)
final_res_Gap_Unif<-matrix(0,12,num_simulation)
for(ii in setting.name){
  wcs = mclapply(1:num_simulation,function(i){
    set.seed(i)
    X<-Sim1(settings = ii)
    x1<-X$x
    gapst_Unif <- K.Clust(x1,method="Gap.Unif",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)$model
    score<-(gapst_Unif$Tab[,3][2]-gapst_Unif$Tab[,3][1])/gapst_Unif$Tab[,4][2]
    #gaprs_Unif <- with(gapst_Unif,maxSE(Tab[,"gap"],Tab[,"SE.sim"]))
    return(score)
  },mc.cores = num_cores)
  final_res_Gap_Unif[which(ii==setting.name),]<-unlist(wcs)
}
#save(final_res_Gap_Unif,file="final_res_Gap_Unif_AUC_extra.Rdata")

final_res_Gap_Unif_Null<-matrix(0,12,num_simulation)
Sim_K_1<-function(dimension,option){
  wcs = mclapply(1:num_simulation,function(i){
    x<-Sim_Null(dimension=dimension,option=option)
    gapst_Unif <- K.Clust(x,method="Gap.Unif",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)$model
    score<-(gapst_Unif$Tab[,3][2]-gapst_Unif$Tab[,3][1])/gapst_Unif$Tab[,4][2]
    return(score)
  },mc.cores = num_cores)
  return(wcs)
}
final_res_Gap_Unif_Null[1,]<-unlist(Sim_K_1(dimension=2,option="Unif"))
final_res_Gap_Unif_Null[2,]<-unlist(Sim_K_1(dimension=3,option="Unif"))
final_res_Gap_Unif_Null[3,]<-unlist(Sim_K_1(dimension=4,option="Unif"))
final_res_Gap_Unif_Null[4,]<-unlist(Sim_K_1(dimension=5,option="Unif"))
final_res_Gap_Unif_Null[5,]<-unlist(Sim_K_1(dimension=10,option="Unif"))
final_res_Gap_Unif_Null[6,]<-unlist(Sim_K_1(dimension=20,option="Unif"))
final_res_Gap_Unif_Null[7,]<-unlist(Sim_K_1(dimension=2,option="Normal"))
final_res_Gap_Unif_Null[8,]<-unlist(Sim_K_1(dimension=3,option="Normal"))
final_res_Gap_Unif_Null[9,]<-unlist(Sim_K_1(dimension=4,option="Normal"))
final_res_Gap_Unif_Null[10,]<-unlist(Sim_K_1(dimension=5,option="Normal"))
final_res_Gap_Unif_Null[11,]<-unlist(Sim_K_1(dimension=10,option="Normal"))
final_res_Gap_Unif_Null[12,]<-unlist(Sim_K_1(dimension=20,option="Normal"))

library(pROC)
AUC_Gap_Unif<-matrix(0,12,12)
rownames(AUC_Gap_Unif)<-c("unif_2","unif_3","unif_4","unif_5","unif_10","unif_20",
                          "normal_2","normal_3","normal_4","normal_5","normal_10","normal_20")
colnames(AUC_Gap_Unif)<-setting.name
#pdf("S4.pdf")
#par(mfrow=c(12,9))
for(i in 1:12){
  for(j in 1:12){
    Gap_1<-as.numeric(final_res_Gap_Unif[j,])
    Gap_null<-as.numeric(final_res_Gap_Unif_Null[i,])
    AUC_Gap_Unif[i,j]<-auc(factor(c(rep("1",num_simulation),rep("0",num_simulation))),c(Gap_1,Gap_null),direction="<")
    #plot.roc(factor(c(rep("1",100),rep("0",100))),c(S4_1,S4_null),direction="<")
  }
  
}
#dev.off()
write.csv(AUC_Gap_Unif,file="AUC_Gap_Unif.csv")

#################################Gap statistic PCA
library(cluster)
final_res_Gap_PCA<-matrix(0,12,num_simulation)
for(ii in setting.name){
  wcs = mclapply(1:num_simulation,function(i){
    set.seed(i)
    X<-Sim1(settings = ii)
    x1<-X$x
    gapst_PCA <- K.Clust(x1,method="Gap.PCA",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)$model
    score<-(gapst_PCA$Tab[,3][2]-gapst_PCA$Tab[,3][1])/gapst_PCA$Tab[,4][2]
    #gaprs_PCA <- with(gapst_PCA,maxSE(Tab[,"gap"],Tab[,"SE.sim"]))
    return(score)
  },mc.cores = num_cores)
  final_res_Gap_PCA[which(ii==setting.name),]<-unlist(wcs)
}
#save(final_res_Gap_PCA,file="final_res_Gap_PCA_AUC_extra.Rdata")

final_res_Gap_PCA_Null<-matrix(0,12,num_simulation)
Sim_K_1<-function(dimension,option){
  wcs = mclapply(1:num_simulation,function(i){
    x<-Sim_Null(dimension=dimension,option=option)
    gapst_PCA <- K.Clust(x,method="Gap.PCA",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)$model
    score<-(gapst_PCA$Tab[,3][2]-gapst_PCA$Tab[,3][1])/gapst_PCA$Tab[,4][2]
    return(score)
  },mc.cores = num_cores)
  return(wcs)
}
final_res_Gap_PCA_Null[1,]<-unlist(Sim_K_1(dimension=2,option="PCA"))
final_res_Gap_PCA_Null[2,]<-unlist(Sim_K_1(dimension=3,option="PCA"))
final_res_Gap_PCA_Null[3,]<-unlist(Sim_K_1(dimension=4,option="PCA"))
final_res_Gap_PCA_Null[4,]<-unlist(Sim_K_1(dimension=5,option="PCA"))
final_res_Gap_PCA_Null[5,]<-unlist(Sim_K_1(dimension=10,option="PCA"))
final_res_Gap_PCA_Null[6,]<-unlist(Sim_K_1(dimension=20,option="PCA"))
final_res_Gap_PCA_Null[7,]<-unlist(Sim_K_1(dimension=2,option="Normal"))
final_res_Gap_PCA_Null[8,]<-unlist(Sim_K_1(dimension=3,option="Normal"))
final_res_Gap_PCA_Null[9,]<-unlist(Sim_K_1(dimension=4,option="Normal"))
final_res_Gap_PCA_Null[10,]<-unlist(Sim_K_1(dimension=5,option="Normal"))
final_res_Gap_PCA_Null[11,]<-unlist(Sim_K_1(dimension=10,option="Normal"))
final_res_Gap_PCA_Null[12,]<-unlist(Sim_K_1(dimension=20,option="Normal"))

library(pROC)
AUC_Gap_PCA<-matrix(0,12,12)
rownames(AUC_Gap_PCA)<-c("PCA_2","PCA_3","PCA_4","PCA_5","PCA_10","PCA_20",
                          "normal_2","normal_3","normal_4","normal_5","normal_10","normal_20")
colnames(AUC_Gap_PCA)<-setting.name
#pdf("S4.pdf")
#par(mfrow=c(12,9))
for(i in 1:12){
  for(j in 1:12){
    Gap_1<-as.numeric(final_res_Gap_PCA[j,])
    Gap_null<-as.numeric(final_res_Gap_PCA_Null[i,])
    AUC_Gap_PCA[i,j]<-auc(factor(c(rep("1",num_simulation),rep("0",num_simulation))),c(Gap_1,Gap_null),direction="<")
    #plot.roc(factor(c(rep("1",100),rep("0",100))),c(S4_1,S4_null),direction="<")
  }
  
}
#dev.off()
write.csv(AUC_Gap_PCA,file="AUC_Gap_PCA.csv")


####################prediction strength
final_res_ps<-matrix(0,12,num_simulation)
for(ii in setting.name){
  wcs = mclapply(1:num_simulation,function(i){
    set.seed(i)
    X<-Sim1(settings = ii)
    x1<-X$x
    res<-K.Clust(x1,method="PS",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)$model
    return(max(res$mean.pred[-1]))
  },mc.cores = num_cores)
  final_res_ps[which(ii==setting.name),]<-unlist(wcs)
}
#save(final_res_ps,file="final_res_ps_AUC_extra.Rdata")

final_res_ps_Null<-matrix(0,12,num_simulation)
Sim_K_1<-function(dimension,option){
  wcs = mclapply(1:num_simulation,function(i){
    x<-Sim_Null(dimension=dimension,option=option)
    res<-K.Clust(x,method="PS",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)$model
    return(max(res$mean.pred[-1]))
  },mc.cores = num_cores)
  return(wcs)
}
final_res_ps_Null[1,]<-unlist(Sim_K_1(dimension=2,option="Unif"))
final_res_ps_Null[2,]<-unlist(Sim_K_1(dimension=3,option="Unif"))
final_res_ps_Null[3,]<-unlist(Sim_K_1(dimension=4,option="Unif"))
final_res_ps_Null[4,]<-unlist(Sim_K_1(dimension=5,option="Unif"))
final_res_ps_Null[5,]<-unlist(Sim_K_1(dimension=10,option="Unif"))
final_res_ps_Null[6,]<-unlist(Sim_K_1(dimension=20,option="Unif"))
final_res_ps_Null[7,]<-unlist(Sim_K_1(dimension=2,option="Normal"))
final_res_ps_Null[8,]<-unlist(Sim_K_1(dimension=3,option="Normal"))
final_res_ps_Null[9,]<-unlist(Sim_K_1(dimension=4,option="Normal"))
final_res_ps_Null[10,]<-unlist(Sim_K_1(dimension=5,option="Normal"))
final_res_ps_Null[11,]<-unlist(Sim_K_1(dimension=10,option="Normal"))
final_res_ps_Null[12,]<-unlist(Sim_K_1(dimension=20,option="Normal"))


library(pROC)
AUC_ps<-matrix(0,12,12)
rownames(AUC_ps)<-c("unif_2","unif_3","unif_4","unif_5","unif_10","unif_20",
                    "normal_2","normal_3","normal_4","normal_5","normal_10","normal_20")
colnames(AUC_ps)<-setting.name
#pdf("ps.pdf")

for(i in 1:12){
  for(j in 1:12){
    ps_1<-as.numeric(final_res_ps[j,])
    ps_null<-as.numeric(final_res_ps_Null[i,])
    AUC_ps[i,j]<-auc(factor(c(rep("1",num_simulation),rep("0",num_simulation))),c(ps_1,ps_null),direction="<")
    #plot.roc(factor(c(rep("1",100),rep("0",100))),c(ps_1,ps_null),direction="<")
  }
  
}
#dev.off()
write.csv(AUC_ps,file="ps_AUC.csv")
#-------------------------------------------------LD
#here we use LD function to output the maximum score. In the package S4, the LD method will not output score just for simplicity.
LD<-function(x, K.try = 2:10,n.resample,cut.off=0.8){
  n = nrow(x)
  sub.n = round(n*0.7)
  result_k<-function(x,K.try){
    result = sapply(K.try, function(K){
      print(K)
      #o.cluster <- kmeans(x, centers = K, nstart = 100)$cluster
      result.subsampling = lapply(1:n.resample, function(b){
        index.subsample = sample(1:nrow(x), sub.n, replace = F)
        #index.subsample = unlist(sapply(1:K,function(i) sample(which(o.cluster == i), 
        #                         round(sum(o.cluster == i)*0.7), replace = F)))
        xb = x[index.subsample,]
        km.out <- kmeans(xb, centers = K, nstart = 100)
        # run sparse k-means
        Cb = km.out$cluster      
        
        group.init = rep(NA, n)
        group.init[index.subsample] = Cb
        consensus.matrix = sapply(1:n, function(i){
          if(i %in% index.subsample){
            as.integer(group.init[i] == group.init)##1==NA is NA
          } else rep(NA, n)
        })
        #consensus.matrix.upper = consensus.matrix[upper.tri(consensus.matrix)]
        return(list(clustering = consensus.matrix))
      })
      
      cluster = sapply(result.subsampling, function(x) x$clustering, simplify = "array")
      
      r.mtx  =  apply(cluster, c(1,2), function(x) mean(x != 0, na.rm = T))###the average concensus matrix of subsample
      
      o.cluster <- kmeans(x, centers = K, nstart = 100)$cluster####kmeans to the whole sample
      
      o.mtx = sapply(1:n, function(i){
        as.integer(o.cluster[i] == o.cluster)
      })###the concensus matrix of whole sample
      
      #r.mtx[which(o.mtx == 0)] <- (1 - r.mtx[which(o.mtx == 0)])
      nr<-nrow(x)
      
      score<-sum(r.mtx[which(o.mtx == 1)])/sum(o.mtx == 1)
      #approximation_specficity<-sum(r.mtx[which(o.mtx == 0)])/sum(o.mtx == 0)
      return(score)
    })
  }
  result_data<-result_k(x,K.try)
  if(length(which(is.na(result_data)))!=0){
    result_data[which(is.na(result_data))]<-0
  }
  maxk<-K.try[which.max(result_data)]
  res<-list(maxk=maxk,score=max(result_data))
  return(res)
}

final_res_LD<-matrix(0,12,num_simulation)
for(ii in setting.name){
  wcs = mclapply(1:num_simulation,function(i){
    set.seed(i)
    X<-Sim1(settings = ii)
    x1<-X$x
    res<-LD(x1, K.try = 2:10,n.resample)
    return(res$score)
  },mc.cores = num_cores)
  final_res_LD[which(ii==setting.name),]<-unlist(wcs)
}
#save(final_res_LD,file="LD_AUC_extra.Rdata")

final_res_LD_Null<-matrix(0,12,num_simulation)
Sim_K_1<-function(dimension,option){
  wcs = mclapply(1:num_simulation,function(i){
    x<-Sim_Null(dimension=dimension,option=option)
    res<-LD(x, K.try = 2:10,n.resample)
    return(res$score)
  },mc.cores = num_cores)
  return(wcs)
}
final_res_LD_Null[1,]<-unlist(Sim_K_1(dimension=2,option="Unif"))
final_res_LD_Null[2,]<-unlist(Sim_K_1(dimension=3,option="Unif"))
final_res_LD_Null[3,]<-unlist(Sim_K_1(dimension=4,option="Unif"))
final_res_LD_Null[4,]<-unlist(Sim_K_1(dimension=5,option="Unif"))
final_res_LD_Null[5,]<-unlist(Sim_K_1(dimension=10,option="Unif"))
final_res_LD_Null[6,]<-unlist(Sim_K_1(dimension=20,option="Unif"))
final_res_LD_Null[7,]<-unlist(Sim_K_1(dimension=2,option="Normal"))
final_res_LD_Null[8,]<-unlist(Sim_K_1(dimension=3,option="Normal"))
final_res_LD_Null[9,]<-unlist(Sim_K_1(dimension=4,option="Normal"))
final_res_LD_Null[10,]<-unlist(Sim_K_1(dimension=5,option="Normal"))
final_res_LD_Null[11,]<-unlist(Sim_K_1(dimension=10,option="Normal"))
final_res_LD_Null[12,]<-unlist(Sim_K_1(dimension=20,option="Normal"))

library(pROC)
AUC_LD<-matrix(0,12,12)
rownames(AUC_LD)<-c("unif_2","unif_3","unif_4","unif_5","unif_10","unif_20",
                    "normal_2","normal_3","normal_4","normal_5","normal_10","normal_20")
colnames(AUC_LD)<-setting.name

for(i in 1:12){
  for(j in 1:12){
    LD_1<-as.numeric(final_res_LD[j,])
    LD_null<-as.numeric(final_res_LD_Null[i,])
    AUC_LD[i,j]<-auc(factor(c(rep("1",num_simulation),rep("0",num_simulation))),c(LD_1,LD_null),direction="<")
    #plot.roc(factor(c(rep("1",100),rep("0",100))),c(LD_1,LD_null),direction="<")
  }
  
}
#dev.off()
write.csv(AUC_LD,file="LD_AUC.csv")

#-------------------------------FW

final_res_FW<-matrix(0,12,num_simulation)
for(ii in setting.name){
  wcs = mclapply(1:num_simulation,function(i){
    set.seed(i)
    X<-Sim1(settings = ii)
    x1<-X$x
    res<-nselectboot(x1,B=n.resample,clustermethod=kmeansCBI,krange=2:10,largeisgood=T)
    return(max(res$stabk[-1]))
  },mc.cores = num_cores)
  final_res_FW[which(ii==setting.name),]<-unlist(wcs)
}


final_res_FW_Null<-matrix(0,12,num_simulation)
Sim_K_1<-function(dimension,option){
  wcs = mclapply(1:num_simulation,function(i){
    x<-Sim_Null(dimension=dimension,option=option)
    res<-nselectboot(x,B=n.resample,clustermethod=kmeansCBI,krange=2:10,largeisgood=T)
    return(max(res$stabk[-1]))
  },mc.cores = num_cores)
  return(wcs)
}
final_res_FW_Null[1,]<-unlist(Sim_K_1(dimension=2,option="Unif"))
final_res_FW_Null[2,]<-unlist(Sim_K_1(dimension=3,option="Unif"))
final_res_FW_Null[3,]<-unlist(Sim_K_1(dimension=4,option="Unif"))
final_res_FW_Null[4,]<-unlist(Sim_K_1(dimension=5,option="Unif"))
final_res_FW_Null[5,]<-unlist(Sim_K_1(dimension=10,option="Unif"))
final_res_FW_Null[6,]<-unlist(Sim_K_1(dimension=20,option="Unif"))
final_res_FW_Null[7,]<-unlist(Sim_K_1(dimension=2,option="Normal"))
final_res_FW_Null[8,]<-unlist(Sim_K_1(dimension=3,option="Normal"))
final_res_FW_Null[9,]<-unlist(Sim_K_1(dimension=4,option="Normal"))
final_res_FW_Null[10,]<-unlist(Sim_K_1(dimension=5,option="Normal"))
final_res_FW_Null[11,]<-unlist(Sim_K_1(dimension=10,option="Normal"))
final_res_FW_Null[12,]<-unlist(Sim_K_1(dimension=20,option="Normal"))


library(pROC)
AUC_FW<-matrix(0,12,12)
rownames(AUC_FW)<-c("unif_2","unif_3","unif_4","unif_5","unif_10","unif_20",
                    "normal_2","normal_3","normal_4","normal_5","normal_10","normal_20")
colnames(AUC_FW)<-setting.name

for(i in 1:12){
  for(j in 1:12){
    FW_1<-as.numeric(final_res_FW[j,])
    FW_null<-as.numeric(final_res_FW_Null[i,])
    AUC_FW[i,j]<-auc(factor(c(rep("1",num_simulation),rep("0",num_simulation))),c(FW_1,FW_null),direction="<")
    #plot.roc(factor(c(rep("1",100),rep("0",100))),c(FW_1,FW_null),direction="<")
  }
  
}
#dev.off()
write.csv(AUC_FW,file="FW_AUC.csv")
