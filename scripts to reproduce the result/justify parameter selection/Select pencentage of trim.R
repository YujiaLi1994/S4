##########################
#Determine number of percentage trim rho
##########################


#Here we use S4.tune function. S4.tune is almost the same as KL.S4.
#but it allows to input multiple time percentage, trying to make the code concise and faster in this script.
#Since rho=5% is already recommended, users are encouraged to use KL.S4 and input rho=5%
#This script is just to show how we conduct simulations to determine best rho
S4.tune<-function(x, K.try = 2:10,n.resample,trim=c(0,0.02,0.05,0.08,0.1,0.15,0.2)){
  n = nrow(x)
  sub.n = round(n*0.5)
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
      # #res_score<-rev(sort(rm1.calc(1:nr)))
      for(i in 1:(nr-1)){
        index.order<-order(rm1.calc(i:nr))
        r.mtx[i:nr,i:nr] = r.mtx[i:nr,i:nr][index.order,
                                            index.order]
        #index[i:nr]<-index[i:nr][index.order]
        o.mtx[i:nr,i:nr] = o.mtx[i:nr,i:nr][index.order,
                                            index.order]
      }
      stat.calc1 = function(i,o.mtx,r.mtx){
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
      
      
      res.K<-rep(0,length(trim))
      for(ind_trim in 1:length(trim)){
        trim1<-trim[ind_trim]
        Num<-round(trim1*nr)
        r.mtx1<-r.mtx[(Num+1):nr,(Num+1):nr]
        o.mtx1<-o.mtx[(Num+1):nr,(Num+1):nr]
        #per<-1-trim
        stat = rev(sapply(1:(dim(o.mtx1)[1]), stat.calc1,o.mtx1,r.mtx1))
        stat<-sort(stat,decreasing = T)
        U.min = mean(stat)
        res.K[ind_trim]<-U.min
      }
      
      return(res.K)
    })
  }
  result_data<-result_k(x,K.try)
  
  return(result_data)
}


num_cores=50
num_simulation=100
n.resample=100

library(MASS)
final_res_S4<-matrix(0,12,7)
numk<-c(3,4,4,4,4,2,2,4,4,4,2,2)
setting.name<-c(paste("A",1:7,sep=""),paste("B",1:5,sep=""))
for(ii in setting.name){
  temp<-matrix(NA,7,num_simulation)
  wcs = mclapply(1:num_simulation,function(i){
    set.seed(i)
    X<-Sim1(settings = ii)
    x1<-X$x
    res<-S4.tune(x1, K.try = 2:10,n.resample,trim=c(0,0.02,0.05,0.08,0.1,0.15,0.2))
    estk<-(2:10)[apply(res,1,function(x){return(max(which(x==max(x))))})]
    return(estk)
  },mc.cores = num_cores)
  for(i in 1:num_simulation){
    temp[,i]<-wcs[[i]]
  }
  final_res_S4[which(ii==setting.name),]<-apply(temp,1,function(x){sum(x==numk[which(ii==setting.name)])})
}

colnames(final_res_S4)<-c(0,0.02,0.05,0.08,0.1,0.15,0.2)
write.csv(final_res_S4,file="S4_Select_Trim.csv")
