score_lambda<-function(lambda,x,k,trim=0.05,n.resample=25){
  print(lambda)
  n = nrow(x)
  sub.n = round(n*0.7)
  result.subsampling = lapply(1:n.resample, function(b){
    index.subsample = sample(1:nrow(x), sub.n, replace = F)
    xb = x[index.subsample,]
    #km.out <- kmeans(xb, centers = K, nstart = 100)
    #hc = hclust(dist(xb), method = "complete", members = NULL)
    #cl<-cmeans(xb,3,m=K)
    #set.seed(b)
    b = KMeansSparseCluster1(xb, K=k, wbounds = lambda, nstart = 20, ###################
                             silent = T, maxiter=6, centers=NULL)
    # run sparse k-means
    #Cb = km.out$cluster      
    #Cb <- cutree(hc, k = K)
    Cb = b[[1]]$Cs
    Fb= b[[1]]$ws
    Fb<-as.numeric(Fb>0)
    group.init = rep(NA, n)
    group.init[index.subsample] = Cb
    consensus.matrix = sapply(1:n, function(i){
      if(i %in% index.subsample){
        as.integer(group.init[i] == group.init)
      } else rep(NA, n)
    })
    #consensus.matrix.upper = consensus.matrix[upper.tri(consensus.matrix)]
    #return(list(clustering = consensus.matrix.upper))
    return(list(clustering = consensus.matrix,feature=Fb))
  })
  cluster = sapply(result.subsampling, function(x) x$clustering, simplify = "array")
  
  r.mtx  =  apply(cluster, c(1,2), function(x) mean(x != 0, na.rm = T))
  #cluster = sapply(result.subsampling, function(x) x$clustering)#matrix, col=25,row=4851
  
  #U = apply(cluster, 1, function(x) mean(x != 0, na.rm = T))
  
  #r.mtx = U
  
  #U.sort = sort(U, decreasing = F); #plot(U.sort)
  
  #U.min1 = -sum(abs(U.sort - 0.5))
  o.result <- KMeansSparseCluster1(x, K=k, wbounds = lambda, nstart = 100, ##############
                                   silent = T, maxiter=6, centers=NULL)
  o.cluster<-o.result[[1]]$Cs
  o.feature<-o.result[[1]]$ws
  o.feature<-as.numeric(o.feature>0)
  consensus.matrix = sapply(1:n, function(i){
    as.integer(o.cluster[i] == o.cluster)
  })
  #o.mtx = as.numeric(consensus.matrix[upper.tri(consensus.matrix)])
  o.mtx = consensus.matrix
  r.mtx[which(o.mtx == 0)] <- (1 - r.mtx[which(o.mtx == 0)])
  nr<-nrow(x)
  rm1.calc = function(index){
    rm1 = vector("numeric")
    r.mtx1 = r.mtx[index,index]
    o.mtx1 = o.mtx[index,index]
    for(i in 1:length(index)){
      if(all(c(1,0)%in%o.mtx1[i,])){
        rm1[i] = (mean(r.mtx1[i, which(o.mtx1[i,] == 1)],na.rm = T) + 
                    mean(r.mtx1[i, which(o.mtx1[i,] == 0)],na.rm = T))-1
      } else if(1%in%o.mtx1[i,]){
        rm1[i] = mean(r.mtx1[i, which(o.mtx1[i,] == 1)],na.rm = T)
      } else {
        rm1[i] = mean(r.mtx1[i, which(o.mtx1[i,] == 0)],na.rm = T)
      }
    }
    return(rm1)
  }
  
  for(i in 1:(nr-1)){
    index.order<-order(rm1.calc(i:nr))
    r.mtx[i:nr,i:nr] = r.mtx[i:nr,i:nr][index.order,
                                        index.order]
    #index[i:nr]<-index[i:nr][index.order]
    o.mtx[i:nr,i:nr] = o.mtx[i:nr,i:nr][index.order,
                                        index.order]
  } 
  
  #-----------------------trim=0
  trim=0
  Num<-round(trim*nr)
  r.mtx<-r.mtx[(Num+1):nr,(Num+1):nr]
  o.mtx<-o.mtx[(Num+1):nr,(Num+1):nr]
  stat.calc1 = function(i){
    nr1<-dim(o.mtx)[1]
    if(all(c(1,0)%in%o.mtx[i,1:nr1])){
      (mean(r.mtx[i,1:nr1][which(o.mtx[i,1:nr1] == 0)]) +
         mean(r.mtx[i,1:nr1][which(o.mtx[i,1:nr1] == 1)])) - 1
    }
    else if(1%in%o.mtx[i,1:nr1]){
      mean(r.mtx[i,1:nr1][which(o.mtx[i,1:nr1] == 1)])
    }
    else {
      mean(r.mtx[i,1:nr1][which(o.mtx[i,1:nr1] == 0)])
    }
  }
  stat = rev(sapply(1:(dim(o.mtx)[1]), stat.calc1))
  stat_youden_0<-mean(sort(stat,decreasing = T))
  
  stat.calc1 = function(i){
    nr1<-dim(o.mtx)[1]
    if(1%in%o.mtx[i,1:nr1]){
      return(mean(r.mtx[i,1:nr1][which(o.mtx[i,1:nr1] == 1)]))
    }else{
      return(NA)
    }
    
  }
  stat = rev(sapply(1:(dim(o.mtx)[1]), stat.calc1))
  stat_sens_0<-mean(sort(stat,decreasing = T))
  
  stat.calc1 = function(i){
    nr1<-dim(o.mtx)[1]
    if(0%in%o.mtx[i,1:nr1]){
      return(mean(r.mtx[i,1:nr1][which(o.mtx[i,1:nr1] == 0)]))
    }else{
      return(NA)
    }
  }
  stat = rev(sapply(1:(dim(o.mtx)[1]), stat.calc1))
  stat_spec_0<-mean(sort(stat,decreasing = T))
  
  #-----------------------trim=0.05
  trim=0.05
  Num<-round(trim*nr)
  r.mtx<-r.mtx[(Num+1):nr,(Num+1):nr]
  o.mtx<-o.mtx[(Num+1):nr,(Num+1):nr]
  stat.calc1 = function(i){
    nr1<-dim(o.mtx)[1]
    if(all(c(1,0)%in%o.mtx[i,1:nr1])){
      (mean(r.mtx[i,1:nr1][which(o.mtx[i,1:nr1] == 0)]) +
         mean(r.mtx[i,1:nr1][which(o.mtx[i,1:nr1] == 1)])) - 1
    }
    else if(1%in%o.mtx[i,1:nr1]){
      mean(r.mtx[i,1:nr1][which(o.mtx[i,1:nr1] == 1)])
    }
    else {
      mean(r.mtx[i,1:nr1][which(o.mtx[i,1:nr1] == 0)])
    }
  }
  stat = rev(sapply(1:(dim(o.mtx)[1]), stat.calc1))
  stat_youden_0.05<-mean(sort(stat,decreasing = T))
  
  stat.calc1 = function(i){
    nr1<-dim(o.mtx)[1]
    if(1%in%o.mtx[i,1:nr1]){
      return(mean(r.mtx[i,1:nr1][which(o.mtx[i,1:nr1] == 1)]))
    }else{
      return(NA)
    }
    
  }
  stat = rev(sapply(1:(dim(o.mtx)[1]), stat.calc1))
  stat_sens_0.05<-mean(sort(stat,decreasing = T))
  
  stat.calc1 = function(i){
    nr1<-dim(o.mtx)[1]
    if(0%in%o.mtx[i,1:nr1]){
      return(mean(r.mtx[i,1:nr1][which(o.mtx[i,1:nr1] == 0)]))
    }else{
      return(NA)
    }
  }
  stat = rev(sapply(1:(dim(o.mtx)[1]), stat.calc1))
  stat_spec_0.05<-mean(sort(stat,decreasing = T))
  
  feature = sapply(result.subsampling, function(x) x$feature, simplify = "matrix")
  
  f.mtx  =  apply(feature, 1, function(x) mean(x != 0, na.rm = T))
  F.mtx = o.feature
  
  f.mtx[which(F.mtx == 0)] <- (1 - f.mtx[which(F.mtx == 0)])
  #index<-order(f.mtx)
  #f.mtx<-sort(f.mtx,decreasing = FALSE)
  F.mtx1<-F.mtx
  
  feature_sentivity<-sum(f.mtx[which(F.mtx1 == 1)])/sum(F.mtx1 == 1)
  if(sum(F.mtx1==0)!=0){
    feature_specficity<-sum(f.mtx[which(F.mtx1 == 0)])/sum(F.mtx1 == 0)
    feature_youden<-feature_sentivity+feature_specficity-1
  }else{
    feature_specficity<-NA
    feature_youden<-feature_sentivity
  }
  
  
  
  #U.min2 = 2 - (mean(r.mtx[which(o.mtx == 1)],na.rm = T) + mean(r.mtx[which(o.mtx == 0)],na.rm =T))
  #U.min3 = 1 - mean(r.mtx[which(o.mtx == 1)],na.rm = T)
  res<-list(stat_youden_0.05=stat_youden_0.05,stat_sens_0.05=stat_sens_0.05,stat_spec_0.05=stat_spec_0.05,
            stat_youden_0=stat_youden_0,stat_sens_0=stat_sens_0,stat_spec_0=stat_spec_0,
            feature_youden=feature_youden,feature_sentivity=feature_sentivity,feature_specficity=feature_specficity)
  return(res)
}







################################
#You need to source(supporting_functions.R) to run this script
################################
num_simulation=50
num_cores=50
n.resample = 100 #resampling times
prop = 0.7 #resample n proportion




#j<-1
q = c(50,200) #number of DE features out of 1000
u = c(0.4,0.6,0.8)

#for (j in 1:5){
#j<-1
S4_Sim_high<-function(j,k){
  wcs<-mclapply(1:num_simulation,function(ii){
    k_vector<-2:7
    h<-1000#DE evidence
    set.seed(ii)
    x<-Sim2(h=h,q=q[j],u=u[k])
    wbounds_list<-list(1)
    for(l in 1:length(k_vector)){
      wbounds_list[[l]] = region.lambda(lam1=3,iteration=40,x,k_vector[l])
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
    result_K_sens_0<-list()
    result_K_spec_0<-list()
    result_K_youden_0<-list()
    result_K_sens_0.05<-list()
    result_K_spec_0.05<-list()
    result_K_youden_0.05<-list()
    result_feature_sens<-list()
    result_feature_spec<-list()
    result_feature_youden<-list()

    for(i1 in 1:length(k_vector)){
      k<-k_vector[i1]
      score1<-list()
      score2<-list()
      score3<-list()
      score4<-list()
      score5<-list()
      score6<-list()
      score7<-list()
      score8<-list()
      score9<-list()
      for(i2 in 1:length(wbounds_list[[i1]])){

        temp<-score_lambda(lambda=wbounds_list[[i1]][i2],x,k,trim=0.05,n.resample=100)
        score1[[i2]]<-temp$stat_sens_0
        score2[[i2]]<-temp$stat_spec_0
        score3[[i2]]<-temp$stat_youden_0
        score4[[i2]]<-temp$stat_sens_0.05
        score5[[i2]]<-temp$stat_spec_0.05
        score6[[i2]]<-temp$stat_youden_0.05
        score7[[i2]]<-temp$feature_sentivity
        score8[[i2]]<-temp$feature_specficity
        score9[[i2]]<-temp$feature_youden
      }
      result_K_sens_0[[i1]]<-score1
      result_K_spec_0[[i1]]<-score2
      result_K_youden_0[[i1]]<-score3
      result_K_sens_0.05[[i1]]<-score4
      result_K_spec_0.05[[i1]]<-score5
      result_K_youden_0.05[[i1]]<-score6
      result_feature_sens[[i1]]<-score7
      result_feature_spec[[i1]]<-score8
      result_feature_youden[[i1]]<-score9
    }
    #
    res<-list(result_K_sens_0=result_K_sens_0,result_K_spec_0=result_K_spec_0,result_K_youden_0=result_K_youden_0,
              result_K_sens_0.05=result_K_sens_0.05,result_K_spec_0.05=result_K_spec_0.05,result_K_youden_0.05=result_K_youden_0.05,
              result_feature_sens=result_feature_sens,result_feature_spec=result_feature_spec,
              result_feature_youden=result_feature_youden,x=x,wbounds_list=wbounds_list)

    #res<-wbounds_list
    #res<-x
    return(res)
    
   
  },mc.cores = num_cores)
  return(wcs)
}
  
result_high<-list()
result_high[[1]]<-S4_Sim_high(1,1)
result_high[[2]]<-S4_Sim_high(1,2)
result_high[[3]]<-S4_Sim_high(1,3)
result_high[[4]]<-S4_Sim_high(2,1)
result_high[[5]]<-S4_Sim_high(2,2)
result_high[[6]]<-S4_Sim_high(2,3)
save(result_high,file="Sim_S4_high_ind.Rdata")
#--------------------------------------about how to combine feature score and clustering score
Jaccard.index = function(TrueFeature, SelectedFeature1){
  u = length(union(which(SelectedFeature1 == 1), which(TrueFeature == 1)))
  is = length(intersect(which(SelectedFeature1 == 1), which(TrueFeature == 1)))
  J1 = is/u
  return(J1)
}
RMSE<-function(x,y){
  return(sqrt(sum((x-y)^2)/length(x)))
}

setwd("~/Desktop/S4_revision/S4 high simulation")
load("Sim_S4_high_ind.Rdata")
#load("Sim_S4_high_ind_wbound_list.Rdata")
#load("Sim_S4_high_ind_x.Rdata")
res.k<-matrix(NA,16,6)
res.rmsek<-matrix(NA,16,6)
res.numfeature<-matrix(NA,16,6)
res.ari<-matrix(NA,16,6)


q<-c(50,50,50,200,200,200)

for(i in 1:6){
  print(paste("i=",i,sep=""))
  score_ari1<-c()
  score_ari2<-c()
  score_ari3<-c()
  score_ari4<-c()
  score_k1<-c()
  score_k2<-c()
  score_k3<-c()
  score_k4<-c()
  score_jaccard1<-c()
  score_jaccard2<-c()
  score_jaccard3<-c()
  score_jaccard4<-c()
  score_num1<-c()
  score_num2<-c()
  score_num3<-c()
  score_num4<-c()
  for(j in 1:50){
    print(paste("j=",j,sep=""))
    x<-result_high[[i]][[j]]$x
    result_k<-list()
    result_lambda<-list()
    temp1<-list()
    temp2<-list()
    temp3<-list()
    temp4<-list()
   wbound_list<-result_high[[i]][[j]]$wbounds_list
    for(k in 1:6){
      result_k[[k]]<-unlist(result_high[[i]][[j]]$result_K_youden_0.05[[k]])
      result_lambda[[k]]<-unlist(result_high[[i]][[j]]$result_feature_youden[[k]])
      #----------------------0.25 clustering and 0.75 feature
      temp1[[k]]<-0.25*result_k[[k]]+0.75*result_lambda[[k]]
      #----------------------0.5 clustering and 0.5 feature
      temp2[[k]]<-0.5*result_k[[k]]+0.5*result_lambda[[k]]
      #----------------------0.75 clustering and 0.25 feature
      temp3[[k]]<-0.75*result_k[[k]]+0.25*result_lambda[[k]]
      #----------------------geometric mean
      temp4[[k]]<-sqrt(result_k[[k]]*result_lambda[[k]])
    }
    score_k<-unlist(lapply(result_k,max))
    est_k<-(2:7)[max(which(score_k==max(score_k)))]
    wbound1<-wbound_list[[est_k-1]][which.max(temp1[[est_k-1]])]
    wbound2<-wbound_list[[est_k-1]][which.max(temp2[[est_k-1]])]
    wbound3<-wbound_list[[est_k-1]][which.max(temp3[[est_k-1]])]
    wbound4<-wbound_list[[est_k-1]][which.max(temp4[[est_k-1]])]
    res1<-KMeansSparseCluster(x,K=est_k,wbounds=wbound1,nstart=50)
    ari1<-adjustedRand(c(rep(1,33),rep(2,33),rep(3,33)),res1[[1]]$Cs)[[2]]
    jaccard1<-Jaccard.index(c(rep(1,q[i]),rep(0,1000-q[i])),as.numeric(res1[[1]]$ws>0))
    num1<-sum(res1[[1]]$ws>0)
    
    res2<-KMeansSparseCluster(x,K=est_k,wbounds=wbound2,nstart=50)
    ari2<-adjustedRand(c(rep(1,33),rep(2,33),rep(3,33)),res2[[1]]$Cs)[[2]]
    jaccard2<-Jaccard.index(c(rep(1,q[i]),rep(0,1000-q[i])),as.numeric(res2[[1]]$ws>0))
    num2<-sum(res2[[1]]$ws>0)
    
    res3<-KMeansSparseCluster(x,K=est_k,wbounds=wbound3,nstart=50)
    ari3<-adjustedRand(c(rep(1,33),rep(2,33),rep(3,33)),res3[[1]]$Cs)[[2]]
    jaccard3<-Jaccard.index(c(rep(1,q[i]),rep(0,1000-q[i])),as.numeric(res3[[1]]$ws>0))
    num3<-sum(res3[[1]]$ws>0)
    
    res4<-KMeansSparseCluster(x,K=est_k,wbounds=wbound4,nstart=50)
    ari4<-adjustedRand(c(rep(1,33),rep(2,33),rep(3,33)),res4[[1]]$Cs)[[2]]
    jaccard4<-Jaccard.index(c(rep(1,q[i]),rep(0,1000-q[i])),as.numeric(res4[[1]]$ws>0))
    num4<-sum(res4[[1]]$ws>0)
    score_ari1[j]<-ari1
    score_ari2[j]<-ari2
    score_ari3[j]<-ari3
    score_ari4[j]<-ari4
    score_k1[j]<-est_k
    score_k2[j]<-est_k
    score_k3[j]<-est_k
    score_k4[j]<-est_k
    score_jaccard1[j]<-jaccard1
    score_jaccard2[j]<-jaccard2
    score_jaccard3[j]<-jaccard3
    score_jaccard4[j]<-jaccard4
    score_num1[j]<-num1
    score_num2[j]<-num2
    score_num3[j]<-num3
    score_num4[j]<-num4
  }
  res.ari[1,i]<-mean(score_ari1)
  res.ari[2,i]<-mean(score_ari2)
  res.ari[3,i]<-mean(score_ari3)
  res.ari[4,i]<-mean(score_ari4)
  res.jaccard[1,i]<-mean(score_jaccard1)
  res.jaccard[2,i]<-mean(score_jaccard2)
  res.jaccard[3,i]<-mean(score_jaccard3)
  res.jaccard[4,i]<-mean(score_jaccard4)
  res.numfeature[1,i]<-mean(num1)
  res.numfeature[2,i]<-mean(num2)
  res.numfeature[3,i]<-mean(num3)
  res.numfeature[4,i]<-mean(num4)
  res.k[1,i]<-RMSE(score_k1,3)
  res.k[2,i]<-RMSE(score_k2,3)
  res.k[3,i]<-RMSE(score_k3,3)
  res.k[4,i]<-RMSE(score_k4,3)
  
  
}
res.ari<-round(res.ari,3)
write.csv(res.ari,file="Combine_two_scores_ind.csv")


#--------------------------------------------------------number of score>0.8
setwd("~/Desktop/S4_revision/S4 high simulation")
load("Sim_S4_high_ind.Rdata")
select_K_1<-rep(NA,6)
for(i in 1:6){
  print(paste("i=",i,sep=""))
  temp<-c()
  for(j in 1:50){
    print(paste("j=",j,sep=""))
    result_k<-list()
    for(k in 1:6){
      result_k[[k]]<-unlist(result_high[[i]][[j]]$result_K_youden_0.05[[k]])
      
    }
    score_k<-max(unlist(lapply(result_k,max)))
    est_k<-(2:7)[max(which(score_k==max(score_k)))]
    if(score_k<0.95){
      temp[j]<-1
    }else{
      temp[j]<-est_k
    }
  
  }
  select_K_1[i]<-sum(temp!=1)
}

#res.ari<-round(select_K_1,3)
#write.csv(select_K_1,file="High_dim_select_1.csv")




#-------------------tune sensitivity and trim
Jaccard.index = function(TrueFeature, SelectedFeature1){
  u = length(union(which(SelectedFeature1 == 1), which(TrueFeature == 1)))
  is = length(intersect(which(SelectedFeature1 == 1), which(TrueFeature == 1)))
  J1 = is/u
  return(J1)
}
RMSE<-function(x,y){
  return(sqrt(sum((x-y)^2)/length(x)))
}

setwd("~/Desktop/S4_revision/S4 high simulation")
load("Sim_S4_high_ind.Rdata")
load("Sim_S4_high_ind_wbound_list.Rdata")
load("Sim_S4_high_ind_x.Rdata")
res.k<-matrix(NA,4,6)
res.jaccard<-matrix(NA,4,6)
res.numfeature<-matrix(NA,4,6)
res.ari<-matrix(NA,4,6)


q<-c(50,50,50,200,200,200)

for(i in 1:6){
  print(paste("i=",i,sep=""))
  score_ari1<-c()
  score_ari2<-c()
  score_ari3<-c()
  score_ari4<-c()
  score_k1<-c()
  score_k2<-c()
  score_k3<-c()
  score_k4<-c()
  score_jaccard1<-c()
  score_jaccard2<-c()
  score_jaccard3<-c()
  score_jaccard4<-c()
  score_num1<-c()
  score_num2<-c()
  score_num3<-c()
  score_num4<-c()
  for(j in 1:50){
    print(paste("j=",j,sep=""))
    x<-result_high[[i]][[j]]$x
    result_k1<-list()
    result_k2<-list()
    result_k3<-list()
    result_k4<-list()
    result_lambda1<-list()
    result_lambda2<-list()
    result_lambda3<-list()
    result_lambda4<-list()
    temp1<-list()
    temp2<-list()
    temp3<-list()
    temp4<-list()
    wbound_list<-result_high[[i]][[j]]$wbounds_list
    for(k in 1:6){
      
      #----------------------no trim and sens
      result_k1[[k]]<-unlist(result_high[[i]][[j]]$result_K_sens_0[[k]])
      result_lambda1[[k]]<-unlist(result_high[[i]][[j]]$result_feature_sens[[k]])
      temp1[[k]]<-0.5*result_k1[[k]]+0.5*result_lambda1[[k]]
      #----------------------trim 5% and sensitivity
      result_k2[[k]]<-unlist(result_high[[i]][[j]]$result_K_sens_0.05[[k]])
      result_lambda2[[k]]<-unlist(result_high[[i]][[j]]$result_feature_sens[[k]])
      temp2[[k]]<-0.5*result_k2[[k]]+0.5*result_lambda2[[k]]
      #----------------------no trim youden
      result_k3[[k]]<-unlist(result_high[[i]][[j]]$result_K_youden_0[[k]])
      result_lambda3[[k]]<-unlist(result_high[[i]][[j]]$result_feature_youden[[k]])
      temp3[[k]]<-0.5*result_k3[[k]]+0.5*result_lambda3[[k]]
      #----------------------trim 5% and youden
      result_k4[[k]]<-unlist(result_high[[i]][[j]]$result_K_youden_0.05[[k]])
      result_lambda4[[k]]<-unlist(result_high[[i]][[j]]$result_feature_youden[[k]])
      temp4[[k]]<-0.5*result_k4[[k]]+0.5*result_lambda4[[k]]
    }
    score_1<-unlist(lapply(result_k1,max))
    score_2<-unlist(lapply(result_k2,max))
    score_3<-unlist(lapply(result_k3,max))
    score_4<-unlist(lapply(result_k4,max))
    est_k1<-(2:7)[max(which(score_1==max(score_1)))]
    est_k2<-(2:7)[max(which(score_2==max(score_2)))]
    est_k3<-(2:7)[max(which(score_3==max(score_3)))]
    est_k4<-(2:7)[max(which(score_4==max(score_4)))]
    wbound1<-wbound_list[[est_k1-1]][which.max(temp1[[est_k1-1]])]
    wbound2<-wbound_list[[est_k2-1]][which.max(temp2[[est_k2-1]])]
    wbound3<-wbound_list[[est_k3-1]][which.max(temp3[[est_k3-1]])]
    wbound4<-wbound_list[[est_k4-1]][which.max(temp4[[est_k4-1]])]
    res1<-KMeansSparseCluster(x,K=est_k1,wbounds=wbound1,nstart=50)
    ari1<-adjustedRand(c(rep(1,33),rep(2,33),rep(3,33)),res1[[1]]$Cs)[[2]]
    jaccard1<-Jaccard.index(c(rep(1,q[i]),rep(0,1000-q[i])),as.numeric(res1[[1]]$ws>0))
    num1<-sum(res1[[1]]$ws>0)
    
    res2<-KMeansSparseCluster(x,K=est_k2,wbounds=wbound2,nstart=50)
    ari2<-adjustedRand(c(rep(1,33),rep(2,33),rep(3,33)),res2[[1]]$Cs)[[2]]
    jaccard2<-Jaccard.index(c(rep(1,q[i]),rep(0,1000-q[i])),as.numeric(res2[[1]]$ws>0))
    num2<-sum(res2[[1]]$ws>0)
    
    res3<-KMeansSparseCluster(x,K=est_k3,wbounds=wbound3,nstart=50)
    ari3<-adjustedRand(c(rep(1,33),rep(2,33),rep(3,33)),res3[[1]]$Cs)[[2]]
    jaccard3<-Jaccard.index(c(rep(1,q[i]),rep(0,1000-q[i])),as.numeric(res3[[1]]$ws>0))
    num3<-sum(res3[[1]]$ws>0)
    
    res4<-KMeansSparseCluster(x,K=est_k4,wbounds=wbound4,nstart=50)
    ari4<-adjustedRand(c(rep(1,33),rep(2,33),rep(3,33)),res4[[1]]$Cs)[[2]]
    jaccard4<-Jaccard.index(c(rep(1,q[i]),rep(0,1000-q[i])),as.numeric(res4[[1]]$ws>0))
    num4<-sum(res4[[1]]$ws>0)
    score_ari1[j]<-ari1
    score_ari2[j]<-ari2
    score_ari3[j]<-ari3
    score_ari4[j]<-ari4
    score_k1[j]<-est_k1
    score_k2[j]<-est_k2
    score_k3[j]<-est_k3
    score_k4[j]<-est_k4
    score_jaccard1[j]<-jaccard1
    score_jaccard2[j]<-jaccard2
    score_jaccard3[j]<-jaccard3
    score_jaccard4[j]<-jaccard4
    score_num1[j]<-num1
    score_num2[j]<-num2
    score_num3[j]<-num3
    score_num4[j]<-num4
  }
  res.ari[1,i]<-mean(score_ari1)
  res.ari[2,i]<-mean(score_ari2)
  res.ari[3,i]<-mean(score_ari3)
  res.ari[4,i]<-mean(score_ari4)
  res.jaccard[1,i]<-mean(score_jaccard1)
  res.jaccard[2,i]<-mean(score_jaccard2)
  res.jaccard[3,i]<-mean(score_jaccard3)
  res.jaccard[4,i]<-mean(score_jaccard4)
  res.numfeature[1,i]<-mean(num1)
  res.numfeature[2,i]<-mean(num2)
  res.numfeature[3,i]<-mean(num3)
  res.numfeature[4,i]<-mean(num4)
  res.k[1,i]<-RMSE(score_k1,3)
  res.k[2,i]<-RMSE(score_k2,3)
  res.k[3,i]<-RMSE(score_k3,3)
  res.k[4,i]<-RMSE(score_k4,3)
  
  
}
#res.ari<-round(res.ari,3)
res<-rbind(res.ari,res.jaccard,res.numfeature,res.k)
res<-round(res,3)
write.csv(res,file="Combine_two_scores_ind_sens_and_trim.csv")

q<-c(50,50,50,200,200,200)
u<-c(0.4,0.6,0.8,0.4,0.6,0.8)
#---------------------------------------------density plot of sensitivity and specificity score

setwd("~/Desktop/S4_revision/S4 high simulation")
par(mfrow=c(2,3))
load("Sim_S4_high_ind.Rdata")
for(i in 1:6){
  res.youden<-rep(NA,50)
  res.sens<-rep(NA,50)
  res.spec<-rep(NA,50)
  for(j in 1:50){
    result<-result_high[[i]][[j]]
    result_k<-list()
    result_sens<-list()
    result_spec<-list()
    for(k in 1:6){
      result_k[[k]]<-unlist(result$result_K_youden_0.05[[k]])
      result_sens[[k]]<-unlist(result$result_K_sens_0.05[[k]])
      result_spec[[k]]<-unlist(result$result_K_spec_0.05[[k]])
    }
    index1<-which.max(unlist(lapply(result_k,max))==max(unlist(lapply(result_k,max))))
    index2<-which.max(result_k[[index1]]==max(result_k[[index1]]))
    res.youden[j]<-result_k[[index1]][index2]
    res.sens[j]<-result_sens[[index1]][index2]
    res.spec[j]<-result_spec[[index1]][index2]
  }
  DF <- data.frame(Sensitivity=res.sens, Specificity=res.spec)
  boxplot(DF, col = rainbow(2, s = 0.5),cex.axis=3,ylim=c(0.4,1.1))
  #axis(side = 1, at = c(2,5), labels = c("sensitivity","specificity"))
  # legend("bottomright", fill = rainbow(2, s = 0.5),
  #        legend = c("Sensitivity","Specificity"),text.font=4,cex=2)
  #plot(density(res.youden),col="black",main="density plot",ylim=c(0,5),xlim=c(0,1))
  #title1<-paste("q=",q[i],", u=",u[i],sep="")
  # plot(density(res.sens),col="black",
  #      ylim=c(0,10),xlim=c(0,1),xlab="Score")
  # lines(density(res.spec),col="red")
  #legend("topleft", legend=c("sensitivity","specificity"), fill=c("black","red"))
}



######################check real data cut.off
setwd("~/Desktop/S4_revision/Detect Null/old real data result")
load("final_result_GSE13159.Rdata")
load("final_result_GSE6891.Rdata")
load("final_result_leaf.Rdata")
load("final_result_su.Rdata")
load("ISOLETs4.Rdata")
load("TCGAPANCANs4.Rdata")
load("SNPs4.Rdata")
load("final_result_Tanbinmat_100resample.Rdata")
