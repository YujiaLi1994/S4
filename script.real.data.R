
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
res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num_cores=10)
result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                  silent = FALSE, maxiter=6, centers=NULL)
#compute ARI
clus_label_gap = result.gap[[1]]$Cs # class label
ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
gap_k<-res.gap$res_k
num_gap<-sum(result.gap[[1]]$ws>0)

set.seed(12315)
res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num_cores=10)
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
res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num_cores=10)
result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                  silent = FALSE, maxiter=6, centers=NULL)
#compute ARI
clus_label_gap = result.gap[[1]]$Cs # class label
ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
gap_k<-res.gap$res_k
num_gap<-sum(result.gap[[1]]$ws>0)

set.seed(12315)
res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num_cores=10)
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
res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num_cores=10)
result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                  silent = FALSE, maxiter=6, centers=NULL)
#compute ARI
clus_label_gap = result.gap[[1]]$Cs # class label
ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
gap_k<-res.gap$res_k
num_gap<-sum(result.gap[[1]]$ws>0)

set.seed(12315)
res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num_cores=10)
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
res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num_cores=10)
result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                  silent = FALSE, maxiter=6, centers=NULL)
#compute ARI
clus_label_gap = result.gap[[1]]$Cs # class label
ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
gap_k<-res.gap$res_k
num_gap<-sum(result.gap[[1]]$ws>0)

set.seed(12315)
res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num_cores=10)
result.ps = KMeansSparseCluster(x, K=res.ps$res_k, wbounds = res.ps$res_w, nstart = 100, ######
                                 silent = FALSE, maxiter=6, centers=NULL)

#compute ARI

clus_label_ps = result.ps[[1]]$Cs # class label
ari_ps<-adjustedRand(true.label,clus_label_ps)[[2]]
ps_k<-res.ps$res_k
num_ps<-sum(result.ps[[1]]$ws>0)####number of feature selected


#end2<-Sys.time()
#end2-start1
load("resleaf.Rdata")
res.s4=final_res$s4
data_K4<-x[,which(result.s4[[1]]$ws>0)]

data_K4<-ds.Leaf
data_K4<-data_K4$data
result.s4 = KMeansSparseCluster(data_K4, K=res.s4$optimal_k, wbounds = res.s4$optimal_lambda, nstart = 100, ######
                                silent = FALSE, maxiter=6, centers=NULL)

data_K4<-data_K4[,which(result.s4[[1]]$ws>0)]
#myBreaks <- seq(-1,1,0.01)
row_distance = dist(data_K4, method = "euclidean")
row_cluster = hclust(row_distance, method = "ward.D")
col_distance = dist(t(data_K4), method = "euclidean")
col_cluster = hclust(col_distance, method = "ward.D")
#par(cex.main=0.8)
#colnames(per_MiRNA)<-cancerType1
heatmap.2(data_K4,
          #RowSideColors = c(    # grouping row-variables into different
          # rep("red", num_sample1),   # categories, Measurement 1-3: green    # Measurement 4-6: blue
          #  rep("blue", num_sample1)),    # Measurement 7-10: red
          main = paste("prolif_pan_analysis",sep=""),
          #xlab=cancerType1,# heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",
          # turns off trace lines inside the heat map
          margins =c(12,9),
          # RowSideColors = c(    # grouping row-variables into different
          #   rep("gray", 12),   # categories, Measurement 1-3: green
          #   rep("blue", 12),    # Measurement 4-6: blue
          #   rep("black", 12)),
          #col=my_palette,# widens margins around plot
          col=greenred,       # use on color palette defined earlier
          dendrogram="both",     # only draw a row dendrogram
          #Rowv=as.dendrogram(row_cluster),
          Colv =as.dendrogram(col_cluster))            # turn off column clustering

##############################


data_K4<-ds.Leaf
label=data_K4$label
data_K4<-data_K4$data
result.s4 = KMeansSparseCluster(data_K4, K=4, wbounds = 10, nstart = 100, ######
                                silent = FALSE, maxiter=6, centers=NULL)
sum(result.s4[[1]]$ws>0)
data_K4<-data_K4[,which(result.s4[[1]]$ws>0)]
#myBreaks <- seq(-1,1,0.01)
row_distance = dist(data_K4, method = "euclidean")
row_cluster = hclust(row_distance, method = "ward.D")
col_distance = dist(t(data_K4), method = "euclidean")
col_cluster = hclust(col_distance, method = "ward.D")
#par(cex.main=0.8)
#colnames(per_MiRNA)<-cancerType1
heatmap.2(data_K4,
          #RowSideColors = c(    # grouping row-variables into different
          # rep("red", num_sample1),   # categories, Measurement 1-3: green    # Measurement 4-6: blue
          #  rep("blue", num_sample1)),    # Measurement 7-10: red
          main = paste("prolif_pan_analysis",sep=""),
          #xlab=cancerType1,# heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",
          # turns off trace lines inside the heat map
          margins =c(12,9),
          RowSideColors = c(rep("gray", 16), rep("blue", 16),rep("black", 16),rep("red",16)),
          col=greenred,       # use on color palette defined earlier
          dendrogram="both",     # only draw a row dendrogram
          #Rowv=as.dendrogram(row_cluster),
          Colv =as.dendrogram(col_cluster))            # turn off column clustering







################################
heatmap.2(data_K4,
          #RowSideColors = c(    # grouping row-variables into different
          # rep("red", num_sample1),   # categories, Measurement 1-3: green    # Measurement 4-6: blue
          #  rep("blue", num_sample1)),    # Measurement 7-10: red
          main = paste("prolif_pan_analysis",sep=""),
          #xlab=cancerType1,# heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",
          # turns off trace lines inside the heat map
          margins =c(12,9),
          RowSideColors = as.character(clus_label_s4),
          #col=my_palette,# widens margins around plot
          col=greenred,       # use on color palette defined earlier
          dendrogram="col",     # only draw a row dendrogram
          #Rowv=as.dendrogram(row_cluster),
          Colv =as.dendrogram(col_cluster))            # turn off column clustering




res<-prcomp(data_K4, center = TRUE,scale. = TRUE)

pca1<-res$x[,1]
pca2<-res$x[,2]
pca3<-res$x[,3]
plot(pca1,pca2)
library("plot3D")
scatter3D(pca1, pca2, pca3, clab = c("Sepal", "Width (cm)"))

save(res.s4,file="k=3GSE47474.Rdata")



d <- dist(data_K4) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=3) # k is the number of dim
fit # view results

# plot solution
x1 <- fit$points[,1]
y1 <- fit$points[,2]
z1 <- fit$points[,3]
scatter3D(x1, y1, z1, clab = c("Sepal", "Width (cm)"))
plot(x1, y1, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS")
text(x1, y1, labels = row.names(mydata), cex=.7)


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
res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num_cores=10)
result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                  silent = FALSE, maxiter=6, centers=NULL)
#compute ARI
clus_label_gap = result.gap[[1]]$Cs # class label
ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
gap_k<-res.gap$res_k
num_gap<-sum(result.gap[[1]]$ws>0)

set.seed(12315)
res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num_cores=10)
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
res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num_cores=10)
result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                  silent = FALSE, maxiter=6, centers=NULL)
#compute ARI
clus_label_gap = result.gap[[1]]$Cs # class label
ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
gap_k<-res.gap$res_k
num_gap<-sum(result.gap[[1]]$ws>0)

set.seed(12315)
res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num_cores=20)
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
res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num_cores=10)
result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                 silent = FALSE, maxiter=6, centers=NULL)
#compute ARI
clus_label_gap = result.gap[[1]]$Cs # class label
ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
gap_k<-res.gap$res_k
num_gap<-sum(result.gap[[1]]$ws>0)

set.seed(12315)
res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num_cores=10)
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
res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num_cores=10)
result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                  silent = FALSE, maxiter=6, centers=NULL)
#compute ARI
clus_label_gap = result.gap[[1]]$Cs # class label
ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
gap_k<-res.gap$res_k
num_gap<-sum(result.gap[[1]]$ws>0)

set.seed(12315)
res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num_cores=10)
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
res.gap<-KL.Gap(x,k_vector=k_vector,lambda_list=wbounds_list,n.perms=n.perms,num_cores=10)
result.gap = KMeansSparseCluster(x, K=res.gap$res_k, wbounds = res.gap$res_w, nstart = 100, ######
                                  silent = FALSE, maxiter=6, centers=NULL)
#compute ARI
clus_label_gap = result.gap[[1]]$Cs # class label
ari_gap<-adjustedRand(true.label,clus_label_gap)[[2]]
gap_k<-res.gap$res_k
num_gap<-sum(result.gap[[1]]$ws>0)

set.seed(12315)
res.ps<-KL.PS(x,lambda_list=wbounds_list,cv=2,k_vector=k_vector,M=n.div,num_cores=20)
result.ps = KMeansSparseCluster(x, K=res.ps$res_k, wbounds = res.ps$res_w, nstart = 100, ######
                                 silent = FALSE, maxiter=6, centers=NULL)

#compute ARI

clus_label_ps = result.ps[[1]]$Cs # class label
ari_ps<-adjustedRand(true.label,clus_label_ps)[[2]]
ps_k<-res.ps$res_k
num_ps<-sum(result.ps[[1]]$ws>0)####number of feature selected

#*********************
final_res<-list(s4=res.s4,ps=res.ps,gap=res.gap) 
save(final_res,file = "res6892.Rdata")
final_res<-list(res.s4=res.s4,res.ps=res.ps,res.gap=res.gap,s4_k=s4_k,ari_s4=ari_s4,num_s4=num_s4,gap_k=gap_k,ari_gap=ari_gap,
                num_gap=num_gap,ps_k=ps_k,ari_ps=ari_ps,num_ps=num_ps)
save(final_res,file = "resISOLET.Rdata")
