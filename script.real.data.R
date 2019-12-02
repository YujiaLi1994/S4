
library(KLClust)
n.resample=500
n.div=500
n.perms=500
#data(ds.GSE13159)
x<-ds.GSE17855$data
#true_lable<-data_tanbin_BRCA$label
#table(true_lable)
true.label<-ds.GSE17855$label
#wbounds_list<-list(1)
h<-dim(x)[2]
k_vector<-c(2,3,4,5,6,7)
set.seed(12315)
start1<-Sys.time()
wbounds_list<-mclapply(1:length(k_vector),function(i){
  K<-k_vector[i]
  error = T
  ind<-0
  while(error){
    result= tryCatch(region.lambda(lam1=3,iteration=30,x,K), error = function(x) NA)
    if(!is.na(result[1])){
      error = F
    }
    ind<-ind+1
    if(ind>3){break}
  }
  return(result)
  
},mc.cores = 20)
#wbounds = region.lamda(b1=6,b2=5,x=x,h=h,K=3)
##test whether it is dense enough
# temp<-KMeansSparseCluster(x,K=2,wbounds=19.3,nstart=100)
# num<-rep(0,length(temp))
# for(i in 1:length(num)){
#   num[i]<-sum(temp[[i]]$ws>0)
# }
for(l in 1:length(k_vector)){
  temp<-KMeansSparseCluster(x,K=k_vector[l],wbounds=wbounds_list[[l]],nstart=100)
  num<-rep(0,length(temp))
  for(i in 1:length(num)){
    num[i]<-sum(temp[[i]]$ws>0)
  }
  if(sum(num==ncol(x))>0){
    wbounds_list[[l]]<-wbounds_list[[l]][1:(min(which(num==ncol(x)))-2)]
  }
}




set.seed(12315)
res.s4 = KL.S4(x, lambda_list=wbounds_list,trim=0.05,k_vector=k_vector,n.resample=n.resample,num.cores=20)

result.s4 = KMeansSparseCluster(x, K=res.s4$optimal_k, wbounds = res.s4$optimal_lambda, nstart = 100, ######
                                silent = FALSE, maxiter=6, centers=NULL)

start2<-Sys.time()
start2-start1

#compute ARI
clus_label_s4 = result.s4[[1]]$Cs # class label
ari_s4<-adjustedRand(true.label,clus_label_s4)[[2]]
s4_k<-res.s4$optimal_k
num_s4<-sum(result.s4[[1]]$ws>0)####numb

set.seed(12315)
res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num.cores=10)
result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                  silent = FALSE, maxiter=6, centers=NULL)
#compute ARI
clus_label_gap = result.gap[[1]]$Cs # class label
ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
gap_k<-res.gap$res_k
num_gap<-sum(result.gap[[1]]$ws>0)

set.seed(12315)
res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num.cores=10)
result.ps = KMeansSparseCluster(x, K=res.ps$res_k, wbounds = res.ps$res_w, nstart = 100, ######
                                 silent = FALSE, maxiter=6, centers=NULL)

#compute ARI

clus_label_ps = result.ps[[1]]$Cs # class label
ari_ps<-adjustedRand(true.label,clus_label_ps)[[2]]
ps_k<-res.ps$res_k
num_ps<-sum(result.ps[[1]]$ws>0)####number of feature selected





###########GSE13159


n.resample=500
n.div=500
n.perms=500
#data(ds.GSE13159)
x<-ds.GSE13159$data
#true_lable<-data_tanbin_BRCA$label
#table(true_lable)
true.label<-ds.GSE13159$label
#wbounds_list<-list(1)
h<-dim(x)[2]
k_vector<-c(2,3,4,5,6,7)
set.seed(12315)
wbounds_list<-mclapply(1:length(k_vector),function(i){
  K<-k_vector[i]
  error = T
  ind<-0
  while(error){
    result= tryCatch(region.lambda(lam1=3,iteration=30,x,K), error = function(x) NA)
    if(!is.na(result[1])){
      error = F
    }
    ind<-ind+1
    if(ind>3){break}
  }
  return(result)
  
},mc.cores = 10)
#wbounds = region.lamda(b1=6,b2=5,x=x,h=h,K=3)
##test whether it is dense enough
# temp<-KMeansSparseCluster(x,K=2,wbounds=19.3,nstart=100)
# num<-rep(0,length(temp))
# for(i in 1:length(num)){
#   num[i]<-sum(temp[[i]]$ws>0)
# }
for(l in 1:length(k_vector)){
  temp<-KMeansSparseCluster(x,K=k_vector[l],wbounds=wbounds_list[[l]],nstart=100)
  num<-rep(0,length(temp))
  for(i in 1:length(num)){
    num[i]<-sum(temp[[i]]$ws>0)
  }
  if(sum(num==ncol(x))>0){
    wbounds_list[[l]]<-wbounds_list[[l]][1:(min(which(num==ncol(x)))-2)]
  }
}
start1<-Sys.time()



set.seed(12315)
res.s4 = KL.S4(x, lambda_list=wbounds_list,trim=0.05,k_vector=k_vector,n.resample=n.resample,num.cores=10)

result.s4 = KMeansSparseCluster(x, K=res.s4$optimal_k, wbounds = res.s4$optimal_lambda, nstart = 100, ######
                                silent = FALSE, maxiter=6, centers=NULL)


#compute ARI
clus_label_s4 = result.s4[[1]]$Cs # class label
ari_s4<-adjustedRand(true.label,clus_label_s4)[[2]]
s4_k<-res.s4$optimal_k
num_s4<-sum(result.s4[[1]]$ws>0)####numb

set.seed(12315)
res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num.cores=10)
result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                  silent = FALSE, maxiter=6, centers=NULL)
#compute ARI
clus_label_gap = result.gap[[1]]$Cs # class label
ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
gap_k<-res.gap$res_k
num_gap<-sum(result.gap[[1]]$ws>0)

set.seed(12315)
res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num.cores=10)
result.ps = KMeansSparseCluster(x, K=res.ps$res_k, wbounds = res.ps$res_w, nstart = 100, ######
                                 silent = FALSE, maxiter=6, centers=NULL)

#compute ARI

clus_label_ps = result.ps[[1]]$Cs # class label
ari_ps<-adjustedRand(true.label,clus_label_ps)[[2]]
ps_k<-res.ps$res_k
num_ps<-sum(result.ps[[1]]$ws>0)####number of feature selected

#########################GSE6891

n.resample=500
n.div=500
n.perms=500
x<-ds.GSE6891$data
#true_lable<-data_tanbin_BRCA$label
#table(true_lable)
true.label<-ds.GSE6891$label
#wbounds_list<-list(1)
h<-dim(x)[2]
k_vector<-c(2,3,4,5,6,7)
set.seed(12315)
wbounds_list<-mclapply(1:length(k_vector),function(i){
  K<-k_vector[i]
  error = T
  ind<-0
  while(error){
    result= tryCatch(region.lambda(lam1=3,iteration=40,x,K), error = function(x) NA)
    if(!is.na(result[1])){
      error = F
    }
    ind<-ind+1
    if(ind>3){break}
  }
  return(result)
  
},mc.cores = 10)
#wbounds = region.lamda(b1=6,b2=5,x=x,h=h,K=3)
##test whether it is dense enough
# temp<-KMeansSparseCluster(x,K=2,wbounds=19.3,nstart=100)
# num<-rep(0,length(temp))
# for(i in 1:length(num)){
#   num[i]<-sum(temp[[i]]$ws>0)
# }
for(l in 1:length(k_vector)){
  temp<-KMeansSparseCluster(x,K=k_vector[l],wbounds=wbounds_list[[l]],nstart=100)
  num<-rep(0,length(temp))
  for(i in 1:length(num)){
    num[i]<-sum(temp[[i]]$ws>0)
  }
  if(sum(num==ncol(x))>0){
    wbounds_list[[l]]<-wbounds_list[[l]][1:(min(which(num==ncol(x)))-2)]
  }
}
start1<-Sys.time()



set.seed(12315)
res.s4 = KL.S4(x, lambda_list=wbounds_list,trim=0.05,k_vector=k_vector,n.resample=n.resample,num.cores=10)

result.s4 = KMeansSparseCluster(x, K=res.s4$optimal_k, wbounds = res.s4$optimal_lambda, nstart = 100, ######
                                silent = FALSE, maxiter=6, centers=NULL)


#compute ARI
clus_label_s4 = result.s4[[1]]$Cs # class label
ari_s4<-adjustedRand(true.label,clus_label_s4)[[2]]
s4_k<-res.s4$optimal_k
num_s4<-sum(result.s4[[1]]$ws>0)####numb

set.seed(12315)
res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num.cores=10)
result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                  silent = FALSE, maxiter=6, centers=NULL)
#compute ARI
clus_label_gap = result.gap[[1]]$Cs # class label
ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
gap_k<-res.gap$res_k
num_gap<-sum(result.gap[[1]]$ws>0)

set.seed(12315)
res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num.cores=10)
result.ps = KMeansSparseCluster(x, K=res.ps$res_k, wbounds = res.ps$res_w, nstart = 100, ######
                                 silent = FALSE, maxiter=6, centers=NULL)

#compute ARI

clus_label_ps = result.ps[[1]]$Cs # class label
ari_ps<-adjustedRand(true.label,clus_label_ps)[[2]]
ps_k<-res.ps$res_k
num_ps<-sum(result.ps[[1]]$ws>0)####number of feature selected





#########################GSE47474

n.resample=500
n.div=500
n.perms=500
x<-ds.GSE47474$data
#true_lable<-data_tanbin_BRCA$label
#table(true_lable)
true.label<-ds.GSE47474$label
#wbounds_list<-list(1)
h<-dim(x)[2]
k_vector<-c(2,3,4,5,6,7)
#k_vector<-c(3)
set.seed(12315)
wbounds_list<-mclapply(1:length(k_vector),function(i){
  K<-k_vector[i]
  error = T
  ind<-0
  while(error){
    result= tryCatch(region.lambda(lam1=3,iteration=40,x,K), error = function(x) NA)
    if(!is.na(result[1])){
      error = F
    }
    ind<-ind+1
    if(ind>3){break}
  }
  return(result)
  
},mc.cores = 10)
#wbounds = region.lamda(b1=6,b2=5,x=x,h=h,K=3)
##test whether it is dense enough
# temp<-KMeansSparseCluster(x,K=2,wbounds=19.3,nstart=100)
# num<-rep(0,length(temp))
# for(i in 1:length(num)){
#   num[i]<-sum(temp[[i]]$ws>0)
# }
for(l in 1:length(k_vector)){
  temp<-KMeansSparseCluster(x,K=k_vector[l],wbounds=wbounds_list[[l]],nstart=100)
  num<-rep(0,length(temp))
  for(i in 1:length(num)){
    num[i]<-sum(temp[[i]]$ws>0)
  }
  if(sum(num==ncol(x))>0){
    wbounds_list[[l]]<-wbounds_list[[l]][1:(min(which(num==ncol(x)))-2)]
  }
}
start1<-Sys.time()



set.seed(12315)
res.s4 = KL.S4(x, lambda_list=wbounds_list,trim=0.05,k_vector=k_vector,n.resample=n.resample,num.cores=10)

#load("k=4GSE47474.Rdata")
result.s4 = KMeansSparseCluster(x, K=res.s4$optimal_k, wbounds = res.s4$optimal_lambda, nstart = 100, ######
                                silent = FALSE, maxiter=6, centers=NULL)


#compute ARI
clus_label_s4 = result.s4[[1]]$Cs # class label
ari_s4<-adjustedRand(true.label,clus_label_s4)[[2]]
s4_k<-res.s4$optimal_k
num_s4<-sum(result.s4[[1]]$ws>0)####numb

set.seed(12315)
res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num.cores=10)
result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                  silent = FALSE, maxiter=6, centers=NULL)
#compute ARI
clus_label_gap = result.gap[[1]]$Cs # class label
ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
gap_k<-res.gap$res_k
num_gap<-sum(result.gap[[1]]$ws>0)

set.seed(12315)
res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num.cores=10)
result.ps = KMeansSparseCluster(x, K=res.ps$res_k, wbounds = res.ps$res_w, nstart = 100, ######
                                 silent = FALSE, maxiter=6, centers=NULL)

#compute ARI

clus_label_ps = result.ps[[1]]$Cs # class label
ari_ps<-adjustedRand(true.label,clus_label_ps)[[2]]
ps_k<-res.ps$res_k
num_ps<-sum(result.ps[[1]]$ws>0)####number of feature selected




#############ds.ISOLET

n.resample=500
n.div=500
n.perms=500
x<-ds.ISOLET$data
#true_lable<-data_tanbin_BRCA$label
#table(true_lable)
true.label<-ds.ISOLET$label
#wbounds_list<-list(1)
h<-dim(x)[2]
k_vector<-c(2,3,4,5,6,7)
set.seed(12315)
wbounds_list<-mclapply(1:length(k_vector),function(i){
  K<-k_vector[i]
  error = T
  ind<-0
  while(error){
    result= tryCatch(region.lambda(lam1=3,iteration=40,x,K), error = function(x) NA)
    if(!is.na(result[1])){
      error = F
    }
    ind<-ind+1
    if(ind>3){break}
  }
  return(result)
  
},mc.cores = 10)
#wbounds = region.lamda(b1=6,b2=5,x=x,h=h,K=3)
##test whether it is dense enough
# temp<-KMeansSparseCluster(x,K=2,wbounds=19.3,nstart=100)
# num<-rep(0,length(temp))
# for(i in 1:length(num)){
#   num[i]<-sum(temp[[i]]$ws>0)
# }
for(l in 1:length(k_vector)){
  temp<-KMeansSparseCluster(x,K=k_vector[l],wbounds=wbounds_list[[l]],nstart=100)
  num<-rep(0,length(temp))
  for(i in 1:length(num)){
    num[i]<-sum(temp[[i]]$ws>0)
  }
  if(sum(num==ncol(x))>0){
    wbounds_list[[l]]<-wbounds_list[[l]][1:(min(which(num==ncol(x)))-2)]
  }
}
start1<-Sys.time()



set.seed(12315)
res.s4 = KL.S4(x, lambda_list=wbounds_list,trim=0.05,k_vector=k_vector,n.resample=n.resample,num.cores=10)

result.s4 = KMeansSparseCluster(x, K=res.s4$optimal_k, wbounds = res.s4$optimal_lambda, nstart = 100, ######
                                silent = FALSE, maxiter=6, centers=NULL)


#compute ARI
clus_label_s4 = result.s4[[1]]$Cs # class label
ari_s4<-adjustedRand(true.label,clus_label_s4)[[2]]
s4_k<-res.s4$optimal_k
num_s4<-sum(result.s4[[1]]$ws>0)####numb

set.seed(12315)
res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num.cores=10)
result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                  silent = FALSE, maxiter=6, centers=NULL)
#compute ARI
clus_label_gap = result.gap[[1]]$Cs # class label
ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
gap_k<-res.gap$res_k
num_gap<-sum(result.gap[[1]]$ws>0)

set.seed(12315)
res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num.cores=10)
result.ps = KMeansSparseCluster(x, K=res.ps$res_k, wbounds = res.ps$res_w, nstart = 100, ######
                                 silent = FALSE, maxiter=6, centers=NULL)

#compute ARI

clus_label_ps = result.ps[[1]]$Cs # class label
ari_ps<-adjustedRand(true.label,clus_label_ps)[[2]]
ps_k<-res.ps$res_k
num_ps<-sum(result.ps[[1]]$ws>0)####number of feature selected








#############ds.pancancer

n.resample=500
n.div=500
n.perms=500
x<-ds.Pancancer$data
#true_lable<-data_tanbin_BRCA$label
#table(true_lable)
true.label<-ds.Pancancer$label

true.label<-true.label[,2]
true.label<-as.character(true.label)
true.label1<-true.label
for(i in 1:length(true.label)){
  if(true.label[i]=="BRCA"){true.label1[i]<-1}
  if(true.label[i]=="COAD"){true.label1[i]<-2}
  if(true.label[i]=="KIRC"){true.label1[i]<-3}
  if(true.label[i]=="LUAD"){true.label1[i]<-4}
  if(true.label[i]=="PRAD"){true.label1[i]<-5}
}
true.label<-true.label1
true.label<-as.numeric(true.label)
#wbounds_list<-list(1)
h<-dim(x)[2]
k_vector<-c(2,3,4,5,6,7)
set.seed(12315)
wbounds_list<-mclapply(1:length(k_vector),function(i){
  K<-k_vector[i]
  error = T
  ind<-0
  while(error){
    result= tryCatch(region.lambda(lam1=3,iteration=40,x,K), error = function(x) NA)
    if(!is.na(result[1])){
      error = F
    }
    ind<-ind+1
    if(ind>3){break}
  }
  return(result)
  
},mc.cores = 20)
#wbounds = region.lamda(b1=6,b2=5,x=x,h=h,K=3)
##test whether it is dense enough
# temp<-KMeansSparseCluster(x,K=2,wbounds=19.3,nstart=100)
# num<-rep(0,length(temp))
# for(i in 1:length(num)){
#   num[i]<-sum(temp[[i]]$ws>0)
# }
for(l in 1:length(k_vector)){
  temp<-KMeansSparseCluster(x,K=k_vector[l],wbounds=wbounds_list[[l]],nstart=100)
  num<-rep(0,length(temp))
  for(i in 1:length(num)){
    num[i]<-sum(temp[[i]]$ws>0)
  }
  if(sum(num==ncol(x))>0){
    wbounds_list[[l]]<-wbounds_list[[l]][1:(min(which(num==ncol(x)))-2)]
  }
}
start1<-Sys.time()



set.seed(12315)
res.s4 = KL.S4(x, lambda_list=wbounds_list,trim=0.05,k_vector=k_vector,n.resample=n.resample,num.cores=20)

result.s4 = KMeansSparseCluster(x, K=res.s4$optimal_k, wbounds = res.s4$optimal_lambda, nstart = 100, ######
                                silent = FALSE, maxiter=6, centers=NULL)


#compute ARI
clus_label_s4 = result.s4[[1]]$Cs # class label
ari_s4<-adjustedRand(true.label,clus_label_s4)[[2]]
s4_k<-res.s4$optimal_k
num_s4<-sum(result.s4[[1]]$ws>0)####numb

set.seed(12315)
res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num.cores=10)
result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                  silent = FALSE, maxiter=6, centers=NULL)
#compute ARI
clus_label_gap = result.gap[[1]]$Cs # class label
ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
gap_k<-res.gap$res_k
num_gap<-sum(result.gap[[1]]$ws>0)

set.seed(12315)
res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num.cores=20)
result.ps = KMeansSparseCluster(x, K=res.ps$res_k, wbounds = res.ps$res_w, nstart = 100, ######
                                 silent = FALSE, maxiter=6, centers=NULL)

#compute ARI

clus_label_ps = result.ps[[1]]$Cs # class label
ari_ps<-adjustedRand(true.label,clus_label_ps)[[2]]
ps_k<-res.ps$res_k
num_ps<-sum(result.ps[[1]]$ws>0)####number of feature selected




#############ds.leaf

n.resample=500
n.div=500
n.perms=500
x<-ds.Leaf$data
#true_lable<-data_tanbin_BRCA$label
#table(true_lable)
true.label<-ds.Leaf$label
#wbounds_list<-list(1)
h<-dim(x)[2]
k_vector<-c(2,3,4,5,6,7)
set.seed(12315)
wbounds_list<-mclapply(1:length(k_vector),function(i){
  K<-k_vector[i]
  error = T
  ind<-0
  while(error){
    result= tryCatch(region.lambda(lam1=3,iteration=40,x,K), error = function(x) NA)
    if(!is.na(result[1])){
      error = F
    }
    ind<-ind+1
    if(ind>3){break}
  }
  return(result)
  
},mc.cores = 10)
#wbounds = region.lamda(b1=6,b2=5,x=x,h=h,K=3)
##test whether it is dense enough
# temp<-KMeansSparseCluster(x,K=2,wbounds=19.3,nstart=100)
# num<-rep(0,length(temp))
# for(i in 1:length(num)){
#   num[i]<-sum(temp[[i]]$ws>0)
# }
for(l in 1:length(k_vector)){
  temp<-KMeansSparseCluster(x,K=k_vector[l],wbounds=wbounds_list[[l]],nstart=100)
  num<-rep(0,length(temp))
  for(i in 1:length(num)){
    num[i]<-sum(temp[[i]]$ws>0)
  }
  if(sum(num==ncol(x))>0){
    wbounds_list[[l]]<-wbounds_list[[l]][1:(min(which(num==ncol(x)))-2)]
  }
}
start1<-Sys.time()



set.seed(12315)
res.s4 = KL.S4(x, lambda_list=wbounds_list,trim=0.05,k_vector=k_vector,n.resample=n.resample,num.cores=10)

result.s4 = KMeansSparseCluster(x, K=res.s4$optimal_k, wbounds = res.s4$optimal_lambda, nstart = 100, ######
                                silent = FALSE, maxiter=6, centers=NULL)


#compute ARI
clus_label_s4 = result.s4[[1]]$Cs # class label
ari_s4<-adjustedRand(true.label,clus_label_s4)[[2]]
s4_k<-res.s4$optimal_k
num_s4<-sum(result.s4[[1]]$ws>0)####numb

set.seed(12315)
res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num.cores=10)
result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                 silent = FALSE, maxiter=6, centers=NULL)
#compute ARI
clus_label_gap = result.gap[[1]]$Cs # class label
ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
gap_k<-res.gap$res_k
num_gap<-sum(result.gap[[1]]$ws>0)

set.seed(12315)
res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num.cores=10)
result.ps = KMeansSparseCluster(x, K=res.ps$res_k, wbounds = res.ps$res_w, nstart = 100, ######
                                silent = FALSE, maxiter=6, centers=NULL)

#compute ARI

clus_label_ps = result.ps[[1]]$Cs # class label
ari_ps<-adjustedRand(true.label,clus_label_ps)[[2]]
ps_k<-res.ps$res_k
num_ps<-sum(result.ps[[1]]$ws>0)####number of feature selected





#############ds.TissueType

n.resample=500
n.div=500
n.perms=500
x<-ds.TissueType$data
#true_lable<-data_tanbin_BRCA$label
#table(true_lable)
true.label<-ds.TissueType$label
#wbounds_list<-list(1)
h<-dim(x)[2]
k_vector<-c(2,3,4,5,6,7)
set.seed(12315)
wbounds_list<-mclapply(1:length(k_vector),function(i){
  K<-k_vector[i]
  error = T
  ind<-0
  while(error){
    result= tryCatch(region.lambda(lam1=3,iteration=40,x,K), error = function(x) NA)
    if(!is.na(result[1])){
      error = F
    }
    ind<-ind+1
    if(ind>3){break}
  }
  return(result)
  
},mc.cores = 10)
#wbounds = region.lamda(b1=6,b2=5,x=x,h=h,K=3)
##test whether it is dense enough
# temp<-KMeansSparseCluster(x,K=2,wbounds=19.3,nstart=100)
# num<-rep(0,length(temp))
# for(i in 1:length(num)){
#   num[i]<-sum(temp[[i]]$ws>0)
# }
for(l in 1:length(k_vector)){
  temp<-KMeansSparseCluster(x,K=k_vector[l],wbounds=wbounds_list[[l]],nstart=100)
  num<-rep(0,length(temp))
  for(i in 1:length(num)){
    num[i]<-sum(temp[[i]]$ws>0)
  }
  if(sum(num==ncol(x))>0){
    wbounds_list[[l]]<-wbounds_list[[l]][1:(min(which(num==ncol(x)))-2)]
  }
}
start1<-Sys.time()



set.seed(12315)
res.s4 = KL.S4(x, lambda_list=wbounds_list,trim=0.05,k_vector=k_vector,n.resample=n.resample,num.cores=10)

result.s4 = KMeansSparseCluster(x, K=res.s4$optimal_k, wbounds = res.s4$optimal_lambda, nstart = 100, ######
                                silent = FALSE, maxiter=6, centers=NULL)


#compute ARI
clus_label_s4 = result.s4[[1]]$Cs # class label
ari_s4<-adjustedRand(true.label,clus_label_s4)[[2]]
s4_k<-res.s4$optimal_k
num_s4<-sum(result.s4[[1]]$ws>0)####numb

set.seed(12315)
res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num.cores=10)
result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                  silent = FALSE, maxiter=6, centers=NULL)
#compute ARI
clus_label_gap = result.gap[[1]]$Cs # class label
ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
gap_k<-res.gap$res_k
num_gap<-sum(result.gap[[1]]$ws>0)

set.seed(12315)
res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num.cores=10)
result.ps = KMeansSparseCluster(x, K=res.ps$res_k, wbounds = res.ps$res_w, nstart = 100, ######
                                 silent = FALSE, maxiter=6, centers=NULL)

#compute ARI

clus_label_ps = result.ps[[1]]$Cs # class label
ari_ps<-adjustedRand(true.label,clus_label_ps)[[2]]
ps_k<-res.ps$res_k
num_ps<-sum(result.ps[[1]]$ws>0)####number of feature selected





#############ds.SNP

n.resample=500
n.div=500
n.perms=500
x<-ds.SNP$data
#true_lable<-data_tanbin_BRCA$label
#table(true_lable)
true.label<-ds.SNP$label
#wbounds_list<-list(1)
h<-dim(x)[2]
k_vector<-c(2,3,4,5,6,7)
set.seed(12315)
wbounds_list<-mclapply(1:length(k_vector),function(i){
  K<-k_vector[i]
  error = T
  ind<-0
  while(error){
    result= tryCatch(region.lambda(lam1=3,iteration=40,x,K), error = function(x) NA)
    if(!is.na(result[1])){
      error = F
    }
    ind<-ind+1
    if(ind>3){break}
  }
  return(result)
  
},mc.cores = 20)
#wbounds = region.lamda(b1=6,b2=5,x=x,h=h,K=3)
##test whether it is dense enough
# temp<-KMeansSparseCluster(x,K=2,wbounds=19.3,nstart=100)
# num<-rep(0,length(temp))
# for(i in 1:length(num)){
#   num[i]<-sum(temp[[i]]$ws>0)
# }
for(l in 1:length(k_vector)){
  temp<-KMeansSparseCluster(x,K=k_vector[l],wbounds=wbounds_list[[l]],nstart=100)
  num<-rep(0,length(temp))
  for(i in 1:length(num)){
    num[i]<-sum(temp[[i]]$ws>0)
  }
  if(sum(num==ncol(x))>0){
    wbounds_list[[l]]<-wbounds_list[[l]][1:(min(which(num==ncol(x)))-2)]
  }
}
start1<-Sys.time()



set.seed(12315)
res.s4 = KL.S4(x, lambda_list=wbounds_list,trim=0.05,k_vector=k_vector,n.resample=n.resample,num.cores=20)

result.s4 = KMeansSparseCluster(x, K=res.s4$optimal_k, wbounds = res.s4$optimal_lambda, nstart = 100, ######
                                silent = FALSE, maxiter=6, centers=NULL)


#compute ARI
clus_label_s4 = result.s4[[1]]$Cs # class label
ari_s4<-adjustedRand(true.label,clus_label_s4)[[2]]
s4_k<-res.s4$optimal_k
num_s4<-sum(result.s4[[1]]$ws>0)####numb

set.seed(12315)
res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num.cores=10)
result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                  silent = FALSE, maxiter=6, centers=NULL)
#compute ARI
clus_label_gap = result.gap[[1]]$Cs # class label
ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
gap_k<-res.gap$res_k
num_gap<-sum(result.gap[[1]]$ws>0)

set.seed(12315)
res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num.cores=20)
result.ps = KMeansSparseCluster(x, K=res.ps$res_k, wbounds = res.ps$res_w, nstart = 100, ######
                                 silent = FALSE, maxiter=6, centers=NULL)

#compute ARI

clus_label_ps = result.ps[[1]]$Cs # class label
ari_ps<-adjustedRand(true.label,clus_label_ps)[[2]]
ps_k<-res.ps$res_k
num_ps<-sum(result.ps[[1]]$ws>0)####number of feature selected


