##' Function to simulate data for 12 settings in simulation I.
##' @title Function to simulate data for 12 settings in simulation I. A1-A7 are well-separated settings and B1-B5 are not-well-separated settings, in terms of ARI and difficulty to separate them from Null data.
##' @param settings Simulate setting from 1 to 10 corresponding to the simulation1 in S4 paper.


##' @return The data matrix corresponding to each setting

##' @export
##' @examples
##' \dontrun{
##' Sim1(settings="2")
##' }



Sim1<-function(settings="A1"){
  if(settings=="A1"){
    x = rbind(matrix(rnorm(50),ncol = 2),
              matrix(c(rnorm(25),(rnorm(25)+5)),ncol = 2),#+5
              matrix(c(rnorm(50)+5,rnorm(50)-3),ncol = 2))#-3
    label=c(rep(1,25),rep(2,25),rep(3,50))
    return(list(x=x,label=label))
  }
  if(settings=="A2"){
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
      x = rbind(m1,m2,m3,m4)
    }
    label=c(rep(1,nrow(m1)),rep(2,nrow(m2)),rep(3,nrow(m3)),rep(4,nrow(m4)))
    return(list(x=x,label=label))
  }

  if(settings=="A3"){
    minm = 0
    while(minm < 1){
      centers = mvrnorm(n = 4, mu = rep(0,5), Sigma = 4*diag(5))
      m1 = matrix(rnorm(5*sample(c(25,50),1)),ncol = 5)
      m1 <- t(apply(m1, 1, function(x) x+centers[1,]))
      m2 = matrix(rnorm(5*sample(c(25,50),1)),ncol = 5)
      m2 <- t(apply(m2, 1, function(x) x+centers[2,]))
      m3 = matrix(rnorm(5*sample(c(25,50),1)),ncol = 5)
      m3 <- t(apply(m3, 1, function(x) x+centers[3,]))
      m4 = matrix(rnorm(5*sample(c(25,50),1)),ncol = 5)
      m4 <- t(apply(m4, 1, function(x) x+centers[4,]))
      m <- as.matrix(dist(rbind(m1,m2,m3,m4)))
      m[1:nrow(m1),1:nrow(m1)] <- 1000
      m[(nrow(m1)+1):(nrow(m1)+nrow(m2)),(nrow(m1)+1):(nrow(m1)+nrow(m2))] <- 1000
      m[(nrow(m1)+nrow(m2)+1):(nrow(m1)+nrow(m2)+nrow(m3)),
        (nrow(m1)+nrow(m2)+1):(nrow(m1)+nrow(m2)+nrow(m3))] <- 1000
      m[(nrow(m1)+nrow(m2)+nrow(m3)+1):(nrow(m1)+nrow(m2)+nrow(m3)+nrow(m4)),
        (nrow(m1)+nrow(m2)+nrow(m3)+1):(nrow(m1)+nrow(m2)+nrow(m3)+nrow(m4))] <- 1000
      minm <- min(m)
      x = rbind(m1,m2,m3,m4)
    }
    label=c(rep(1,nrow(m1)),rep(2,nrow(m2)),rep(3,nrow(m3)),rep(4,nrow(m4)))
    return(list(x=x,label=label))
  }

  if(settings=="A4"){
    minm = 0
    while(minm < 1){
      centers = mvrnorm(n = 4, mu = rep(0,8), Sigma = 3*diag(8))
      m1 = matrix(rnorm(8*sample(c(25,50),1)),ncol = 8)
      m1 <- t(apply(m1, 1, function(x) x+centers[1,]))
      m2 = matrix(rnorm(8*sample(c(25,50),1)),ncol = 8)
      m2 <- t(apply(m2, 1, function(x) x+centers[2,]))
      m3 = matrix(rnorm(8*sample(c(25,50),1)),ncol = 8)
      m3 <- t(apply(m3, 1, function(x) x+centers[3,]))
      m4 = matrix(rnorm(8*sample(c(25,50),1)),ncol = 8)
      m4 <- t(apply(m4, 1, function(x) x+centers[4,]))
      m <- as.matrix(dist(rbind(m1,m2,m3,m4)))
      m[1:nrow(m1),1:nrow(m1)] <- 1000
      m[(nrow(m1)+1):(nrow(m1)+nrow(m2)),(nrow(m1)+1):(nrow(m1)+nrow(m2))] <- 1000
      m[(nrow(m1)+nrow(m2)+1):(nrow(m1)+nrow(m2)+nrow(m3)),
        (nrow(m1)+nrow(m2)+1):(nrow(m1)+nrow(m2)+nrow(m3))] <- 1000
      m[(nrow(m1)+nrow(m2)+nrow(m3)+1):(nrow(m1)+nrow(m2)+nrow(m3)+nrow(m4)),
        (nrow(m1)+nrow(m2)+nrow(m3)+1):(nrow(m1)+nrow(m2)+nrow(m3)+nrow(m4))] <- 1000
      minm <- min(m)
      x = rbind(m1,m2,m3,m4)
    }
    label=c(rep(1,nrow(m1)),rep(2,nrow(m2)),rep(3,nrow(m3)),rep(4,nrow(m4)))
    return(list(x=x,label=label))
  }




  if(settings=="A5"){
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
      x = rbind(m1,m2,m3,m4)
    }
    label=c(rep(1,nrow(m1)),rep(2,nrow(m2)),rep(3,nrow(m3)),rep(4,nrow(m4)))
    return(list(x=x,label=label))
  }
  if(settings=="A6"){
    x = rbind(matrix(c(seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1),
                       seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1),
                       seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)),ncol = 3),
              matrix(c(seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)+10,
                       seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)+10,
                       seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)+10),ncol = 3))
    label=c(rep(1,100),rep(2,100))
    return(list(x=x,label=label))
  }
  if(settings=="A7"){
    x = rbind(matrix(c(seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1),
                       seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1),
                       seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)),ncol = 3),
              matrix(c(seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)+1,
                       seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)+1,
                       seq(-0.5,0.5,length.out = 100)+rnorm(100,sd = 0.1)+1),ncol = 3))
    label=c(rep(1,100),rep(2,100))
    return(list(x=x,label=label))
  }


  if(settings=="B1"){
    x = rbind(matrix(rnorm(50),ncol = 2),
              matrix(c(rnorm(25),rnorm(25)+2.5),ncol = 2),
              matrix(c(rnorm(25)+2.5,rnorm(25)),ncol = 2),
              matrix(c(rnorm(25)+2.5,rnorm(25)+2.5),ncol = 2))
    label=c(rep(1,25),rep(2,25),rep(3,25),rep(4,25))
    return(list(x=x,label=label))
  }
  if(settings=="B2"){
    x = rbind(matrix(rnorm(50),ncol = 2),
              matrix(c(rnorm(25),rnorm(25)+3),ncol = 2),
              matrix(c(rnorm(25)+3,rnorm(25)),ncol = 2),
              matrix(c(rnorm(25)+3,rnorm(25)+3),ncol = 2))
    label=c(rep(1,25),rep(2,25),rep(3,25),rep(4,25))
    return(list(x=x,label=label))
  }
  if(settings=="B3"){
    x = rbind(matrix(rnorm(50),ncol = 2),
              matrix(c(rnorm(25),rnorm(25)+3.5),ncol = 2),
              matrix(c(rnorm(25)+3.5,rnorm(25)),ncol = 2),
              matrix(c(rnorm(25)+3.5,rnorm(25)+3.5),ncol = 2))
    label=c(rep(1,25),rep(2,25),rep(3,25),rep(4,25))
    return(list(x=x,label=label))
  }

  if(settings=="B4"){

    shift<-2
    centers1<-rep(0,5)
    centers2<-c(rep(shift,1),rep(0,4))
    m1<-mvrnorm(n=50,centers1,Sigma = diag(5))
    m2<-mvrnorm(n=50,centers2,Sigma = diag(5))
    x = rbind(m1,m2)
    label=c(rep(1,nrow(m1)),rep(2,nrow(m2)))
    return(list(x=x,label=label))
  }


  if(settings=="B5"){
    shift<-2
    centers1<-rep(0,10)
    centers2<-c(rep(shift,1),rep(0,9))
    m1<-mvrnorm(n=50,centers1,Sigma = diag(10))
    m2<-mvrnorm(n=50,centers2,Sigma = diag(10))
    x = rbind(m1,m2)
    label=c(rep(1,nrow(m1)),rep(2,nrow(m2)))
    return(list(x=x,label=label))
  }
}
