library(S4)
#############################
#Estimate K for Kmeans
##############################
#simulate the data of simulation I by Sim1 function.
data<-Sim1(settings = "A1")#A1-A7 are well separated data and B1-B7 are not-well separated data
#using S4
#trim.S4 is the percentage of trimmed scattered points.
#Cutoff is the threshold to detect the K=1
#n.resample is number of subsampling
res.S4<-K.Clust(data$x,method="S4",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=50)
#users can also implement other methods, see the document and reference for details.
#Take prediction strength for an example
res.PS<-K.Clust(data$x,method="PS",Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=50)

#############################
#Estimate K and lambda simultaneously for sparse Kmeans
##############################
#generate simulation II data
data<-Sim2(h=200,q = 50,u=0.8)#For demo purpose, we use 200 features in total for fast result. 
#using Efficient Algorithm for Choosing Grids of lambda for each K
k_vector<-2:7#search K from 2 to 7
wbounds_list<-list(1)
for(l in 1:length(k_vector)){#for each K, using the algorithm to get 20 lambda.
  wbounds_list[[l]] = region.lambda(lam1=1.5,iteration=20,data,k_vector[l])
}

##This part of code is to delete the lambda which selecting all the genes.
#since our S4 method calculate specificity, the lambda selecting all the genes must be removed.
for(l in 1:length(k_vector)){
  temp<-KMeansSparseCluster(data,K=k_vector[l],wbounds=wbounds_list[[l]],nstart=100)
  num<-rep(0,length(temp))
  for(i in 1:length(num)){
    num[i]<-sum(temp[[i]]$ws>0)
  }
  if(sum(num==ncol(data))>0){
    wbounds_list[[l]]<-wbounds_list[[l]][1:(min(which(num==ncol(data)))-3)]
  }
}

#run S4 method
res.S4<-KL.S4(x=data,lambda_list = wbounds_list,k_vector = k_vector,trim =0.05,n.resample = 50,num.cores = 1)
#run extended Gap statistic method
res.Gap<-KL.Gap(x=data,lambda_list = wbounds_list,k_vector = k_vector,n.perm = 50,num.cores = 1)
#run extended Prediction strength method
res.PS<-KL.PS(x=data,lambda_list = wbounds_list,k_vector = k_vector,cv = 2,M=20,num.cores = 1,cutoff = 0.8)

#############################
#access real data
##############################
data()
#Data sets in package ‘S4’:
# ds.GSE13159                  Microarray data of Leukemia after preprocessing(GSE13159)
# ds.GSE17855                  Microarray data of leukemia after preprocessing(GSE17855)
# ds.GSE47474                  RNA Sequencing data of rat brain after
# preprocessing(GSE47474)
# ds.GSE6891                   Microarray data of leukemia after preprocessing(GSE6891)
# ds.ISOLET                    letter-name dataset after preprocessing
# ds.Leaf                      Plant species leaves dataset after preprocessing
# ds.Pancancer                 RNA Sequencing data after preprocessing(Pancancer)
# ds.SNP                       SNP data after preprocessing
# ds.TissueType                Microarray data after preprocessing(Mammalian tissue
#                                                                  types dataset:)
data<-ds.GSE47474$data
label<-ds.GSE47474$label

