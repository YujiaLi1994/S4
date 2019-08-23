library(KLClust)
Number_of_Simulation = 100
K.try = 2:10 #K to try
n.resample = 100#resampling times
prop = 0.7 #resample n proportion
core=30

#---------------
#simulate the data,code for three summary statistics method
#---------------











######################### sensitivity analysis
wc<-mclapply(1:Number_of_Simulation,function(pp){
  x1 = matrix(runif(2000),ncol = 10)
  x2 = rbind(matrix(rnorm(50),ncol = 2),
             matrix(c(rnorm(25),(rnorm(25)+5)),ncol = 2),#+5
             matrix(c(rnorm(50)+5,rnorm(50)-3),ncol = 2))#-3
  minm = 0
  while(minm < 1){
    centers = mvrnorm(n = 4, mu = rep(0,3), Sigma = 5*diag(3))
    m1 = matrix(rnorm(3*sample(c(25,50),1)),ncol = 3)
    m1 <- t(apply(m1, 1, function(x) x+centers[1,]))
    m2 = matrix(rnorm(3*sample(c(25,50),1)),ncol = 3)
    m2 <- t(apply(m2, 1, function(x) x+centers[2,]))
    m3 = matrix(rnorm(3*sample(c(25,50),1)),ncol = 3)
    m3 <- t(apply(m3, 1, function(x) x+centers[3,]))
    m4 = matrix(rnorm(3*sample(c(25,50),1)),ncol = 3)
    m4 <- t(apply(m4, 1, function(x) x+centers[4,]))
    m <- as.matrix(dist(rbind(m1,m2,m3,m4)))
    m[1:nrow(m1),1:nrow(m1)] <- 1000
    m[(nrow(m1)+1):(nrow(m1)+nrow(m2)),(nrow(m1)+1):(nrow(m1)+nrow(m2))] <- 1000
    m[(nrow(m1)+nrow(m2)+1):(nrow(m1)+nrow(m2)+nrow(m3)),
      (nrow(m1)+nrow(m2)+1):(nrow(m1)+nrow(m2)+nrow(m3))] <- 1000
    m[(nrow(m1)+nrow(m2)+nrow(m3)+1):(nrow(m1)+nrow(m2)+nrow(m3)+nrow(m4)),
      (nrow(m1)+nrow(m2)+nrow(m3)+1):(nrow(m1)+nrow(m2)+nrow(m3)+nrow(m4))] <- 1000
    minm <- min(m)
    x3 = rbind(m1,m2,m3,m4)
  }
  minm = 0
  while(minm < 1){
    centers = mvrnorm(n = 4, mu = rep(0,10), Sigma = 1.9*diag(10))
    m1 = matrix(rnorm(10*sample(c(25,50),1)),ncol = 10)
    m1 <- t(apply(m1, 1, function(x) x+centers[1,]))
    m2 = matrix(rnorm(10*sample(c(25,50),1)),ncol = 10)
    m2 <- t(apply(m2, 1, function(x) x+centers[2,]))
    m3 = matrix(rnorm(10*sample(c(25,50),1)),ncol = 10)
    m3 <- t(apply(m3, 1, function(x) x+centers[3,]))
    m4 = matrix(rnorm(10*sample(c(25,50),1)),ncol = 10)
    m4 <- t(apply(m4, 1, function(x) x+centers[4,]))
    m <- as.matrix(dist(rbind(m1,m2,m3,m4)))
    m[1:nrow(m1),1:nrow(m1)] <- 1000
    m[(nrow(m1)+1):(nrow(m1)+nrow(m2)),(nrow(m1)+1):(nrow(m1)+nrow(m2))] <- 1000
    m[(nrow(m1)+nrow(m2)+1):(nrow(m1)+nrow(m2)+nrow(m3)),
      (nrow(m1)+nrow(m2)+1):(nrow(m1)+nrow(m2)+nrow(m3))] <- 1000
    m[(nrow(m1)+nrow(m2)+nrow(m3)+1):(nrow(m1)+nrow(m2)+nrow(m3)+nrow(m4)),
      (nrow(m1)+nrow(m2)+nrow(m3)+1):(nrow(m1)+nrow(m2)+nrow(m3)+nrow(m4))] <- 1000
    minm <- min(m)
    x4 = rbind(m1,m2,m3,m4)
  }
  x5 = rbind(matrix(c(seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1),
                      seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1),
                      seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)),ncol = 3),
             matrix(c(seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)+10,
                      seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)+10,
                      seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)+10),ncol = 3))
  x6 = rbind(matrix(rnorm(50),ncol = 2),
             matrix(c(rnorm(25),rnorm(25)+2.5),ncol = 2),
             matrix(c(rnorm(25)+2.5,rnorm(25)),ncol = 2),
             matrix(c(rnorm(25)+2.5,rnorm(25)+2.5),ncol = 2))
  x7 = rbind(matrix(rnorm(50),ncol = 2),
             matrix(c(rnorm(25),rnorm(25)+3),ncol = 2),
             matrix(c(rnorm(25)+3,rnorm(25)),ncol = 2),
             matrix(c(rnorm(25)+3,rnorm(25)+3),ncol = 2))
  x8 = rbind(matrix(rnorm(50),ncol = 2),
             matrix(c(rnorm(25),rnorm(25)+3.5),ncol = 2),
             matrix(c(rnorm(25)+3.5,rnorm(25)),ncol = 2),
             matrix(c(rnorm(25)+3.5,rnorm(25)+3.5),ncol = 2))
  
  
  x9 = rbind(matrix(c(seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1),
                      seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1),
                      seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)),ncol = 3),
             matrix(c(seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)+1,
                      seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)+1,
                      seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)+1),ncol = 3))
  x10 = rbind(matrix(c(seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1),
                       seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1),
                       seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)),ncol = 3),
              matrix(c(seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)+1,
                       seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1),
                       seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)),ncol = 3))
  trim<-0
  wc1<-K.Clust(x1, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc2<-K.Clust(x2, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc3<-K.Clust(x3, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc4<-K.Clust(x4, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc5<-K.Clust(x5, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc6<-K.Clust(x6, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc7<-K.Clust(x7, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc8<-K.Clust(x8, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc9<-K.Clust(x9, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc10<-K.Clust(x10, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  awc1<-list(wc1=wc1,wc2=wc2,wc3=wc3,wc4=wc4,wc5=wc5,wc6=wc6,wc7=wc7,wc8=wc8,wc9=wc9,wc10=wc10)
  trim<-0.02
  wc1<-K.Clust(x1, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc2<-K.Clust(x2, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc3<-K.Clust(x3, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc4<-K.Clust(x4, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc5<-K.Clust(x5, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc6<-K.Clust(x6, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc7<-K.Clust(x7, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc8<-K.Clust(x8, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc9<-K.Clust(x9, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc10<-K.Clust(x10, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  
  awc2<-list(wc1=wc1,wc2=wc2,wc3=wc3,wc4=wc4,wc5=wc5,wc6=wc6,wc7=wc7,wc8=wc8,wc9=wc9,wc10=wc10)
  
  trim<-0.05
  wc1<-K.Clust(x1, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc2<-K.Clust(x2, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc3<-K.Clust(x3, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc4<-K.Clust(x4, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc5<-K.Clust(x5, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc6<-K.Clust(x6, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc7<-K.Clust(x7, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc8<-K.Clust(x8, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc9<-K.Clust(x9, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc10<-K.Clust(x10, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  awc3<-list(wc1=wc1,wc2=wc2,wc3=wc3,wc4=wc4,wc5=wc5,wc6=wc6,wc7=wc7,wc8=wc8,wc9=wc9,wc10=wc10)
  
  trim<-0.08
  wc1<-K.Clust(x1, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc2<-K.Clust(x2, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc3<-K.Clust(x3, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc4<-K.Clust(x4, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc5<-K.Clust(x5, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc6<-K.Clust(x6, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc7<-K.Clust(x7, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc8<-K.Clust(x8, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc9<-K.Clust(x9, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc10<-K.Clust(x10, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  awc4<-list(wc1=wc1,wc2=wc2,wc3=wc3,wc4=wc4,wc5=wc5,wc6=wc6,wc7=wc7,wc8=wc8,wc9=wc9,wc10=wc10)
  
  trim<-0.10
  wc1<-K.Clust(x1,Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc2<-K.Clust(x2, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc3<-K.Clust(x3, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc4<-K.Clust(x4, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc5<-K.Clust(x5, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc6<-K.Clust(x6, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc7<-K.Clust(x7, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc8<-K.Clust(x8, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc9<-K.Clust(x9, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc10<-K.Clust(x10, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  awc5<-list(wc1=wc1,wc2=wc2,wc3=wc3,wc4=wc4,wc5=wc5,wc6=wc6,wc7=wc7,wc8=wc8,wc9=wc9,wc10=wc10)
  
  trim<-0.15
  wc1<-K.Clust(x1, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc2<-K.Clust(x2, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc3<-K.Clust(x3, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc4<-K.Clust(x4, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc5<-K.Clust(x5, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc6<-K.Clust(x6, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc7<-K.Clust(x7, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc8<-K.Clust(x8, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc9<-K.Clust(x9, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc10<-K.Clust(x10, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  awc6<-list(wc1=wc1,wc2=wc2,wc3=wc3,wc4=wc4,wc5=wc5,wc6=wc6,wc7=wc7,wc8=wc8,wc9=wc9,wc10=wc10)
  
  trim<-0.20
  wc1<-K.Clust(x1, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc2<-K.Clust(x2, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc3<-K.Clust(x3, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc4<-K.Clust(x4, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc5<-K.Clust(x5, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc6<-K.Clust(x6, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc7<-K.Clust(x7, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc8<-K.Clust(x8, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc9<-K.Clust(x9, Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  wc10<-K.Clust(x10,Kmin=2,Kmax=10,trim.S4=0.05,cutoff=0.8,n.resample=n.resample)
  awc7<-list(wc1=wc1,wc2=wc2,wc3=wc3,wc4=wc4,wc5=wc5,wc6=wc6,wc7=wc7,wc8=wc8,wc9=wc9,wc10=wc10)
  return(list(awc1,awc2,awc3,awc4,awc5,awc6,awc7))
},mc.cores = core)
res1<-rep(0,100)
res2<-rep(0,100)
res3<-rep(0,100)
res4<-rep(0,100)
res5<-rep(0,100)
res6<-rep(0,100)
res7<-rep(0,100)
res8<-rep(0,100)
res9<-rep(0,100)
res10<-rep(0,100)

#res12<-rep(0,100)
result<-matrix(rep(0,10*7),ncol=10)
for(j in 1:7){
  for(i in 1:100){
    res1[i]<-wc[[i]][[j]][[1]]
    res2[i]<-wc[[i]][[j]][[2]]
    res3[i]<-wc[[i]][[j]][[3]]
    res4[i]<-wc[[i]][[j]][[4]]
    res5[i]<-wc[[i]][[j]][[5]]
    res6[i]<-wc[[i]][[j]][[6]]
    res7[i]<-wc[[i]][[j]][[7]]
    res8[i]<-wc[[i]][[j]][[8]]
    res9[i]<-wc[[i]][[j]][[9]]
    res10[i]<-wc[[i]][[j]][[10]]
    
  }
  res1<-sum(res1==1)
  res2<-sum(res2==3)
  res3<-sum(res3==4)
  res4<-sum(res4==4)
  res5<-sum(res5==2)
  res6<-sum(res6==4)
  res7<-sum(res7==4)
  res8<-sum(res8==4)
  res9<-sum(res9==2)
  res10<-sum(res10==2)
  
  result[j,]<-cbind(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10)
}
