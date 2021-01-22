
weighted_l2<-function(x,y,w){
  res<-0
  for(i in 1:length(x)){
    res<-res+w[i]*(x[i]-y[i])^2
  }
  return(sqrt(res))
}






KmeansSparseCluster.permute1<-function(x, K=NULL,  nperms=25, wbounds=NULL,silent=FALSE, nvals=10, centers=NULL){
  if(is.null(wbounds)) wbounds <- exp(seq(log(1.2), log(sqrt(ncol(x))*.9), len=nvals))
  if(min(wbounds) <= 1) stop("Wbounds should be greater than 1, since otherwise only one weight will be nonzero.")
  if(length(wbounds)<2) stop("Wbounds should be a vector of at least two elements.")
  # was seq(1.2, sqrt(ncol(x))*.6, len=10)
  if(is.null(K) && is.null(centers)) stop("Must provide either K or centers.")
  if(!is.null(K) && !is.null(centers)){
    if(nrow(centers)!=K) stop("If K and centers both are provided, then nrow(centers) must equal K!!!")
    if(nrow(centers)==K) K <- NULL
  }
  if(!is.null(centers) && ncol(centers)!=ncol(x)) stop("If centers is provided, then ncol(centers) must equal ncol(x).")
  permx <- list()
  nnonzerows <- NULL
  for(i in 1:nperms){
    permx[[i]] <- matrix(NA, nrow=nrow(x), ncol=ncol(x))
    for(j in 1:ncol(x)) permx[[i]][,j] <- sample(x[,j])
  }
  tots <- NULL
  out <- KMeansSparseCluster1(x, K, wbounds=wbounds, silent=silent, centers=centers)
  for(i in 1:length(out)){
    nnonzerows <- c(nnonzerows, sum(out[[i]]$ws!=0))
    bcss <- GetWCSS(x,out[[i]]$Cs)$bcss.perfeature
    tots <- c(tots, sum(out[[i]]$ws*bcss))
  }
  permtots <- matrix(NA, nrow=length(wbounds), ncol=nperms)
  for(k in 1:nperms){
    if(!silent) cat("Permutation ", k, "of ", nperms, fill=TRUE)
    perm.out <- KMeansSparseCluster1(permx[[k]], K, wbounds=wbounds, silent=silent, centers=centers)
    for(i in 1:length(perm.out)){
      perm.bcss <- GetWCSS(permx[[k]], perm.out[[i]]$Cs)$bcss.perfeature
      permtots[i,k] <- sum(perm.out[[i]]$ws*perm.bcss)
    }
  }
  gaps <- (log(tots)-apply(log(permtots),1,mean))
  out <- list(tots=tots, permtots=permtots, nnonzerows=nnonzerows, gaps=gaps, sdgaps=apply(log(permtots),1,sd), wbounds=wbounds, bestw=wbounds[which.max(gaps)])
  if(!silent) cat(fill=TRUE)
  class(out) <- "KmeansSparseCluster.permute"
  return(out)
}

print.KmeansSparseCluster.permute <- function(x,...){
  cat("Tuning parameter selection results for Sparse K-means Clustering:", fill=TRUE)
  mat <- round(cbind(x$wbounds, x$nnonzerows, x$gaps, x$sdgaps),4)
  dimnames(mat) <- list(1:length(x$wbounds), c("Wbound", "# Non-Zero W's", "Gap Statistic", "Standard Deviation"))
  print(mat, quote=FALSE)
  cat("Tuning parameter that leads to largest Gap statistic: ", x$bestw, fill=TRUE)
}

plot.KmeansSparseCluster.permute <- function(x,...){
  plot(x$nnonzerows, x$gaps, log="x", main="Gap Statistics", xlab="# Non-zero Wj's", ylab="")
  lines(x$nnonzerows, x$gaps)
}

KMeansSparseCluster1<-function(x, K=NULL, wbounds=NULL, nstart=20, silent=FALSE, maxiter=6, centers=NULL){
  # The criterion is : minimize_{w, C} sum_j w_j (WCSS_j - TSS_j) s.t. ||w||_2=1, ||w||_1<=s, w_j>=0
  # x is the data, nxp
  # K is the number of clusters desired
  # wbounds is a vector of L1 constraints on w, of the form  sum(abs(w))<=wbounds[i]
  if(is.null(K) && is.null(centers)) stop("Must provide either K or centers.")
  if(!is.null(K) && !is.null(centers)){
    if(nrow(centers)!=K) stop("If K and centers both are provided, then nrow(centers) must equal K!!!")
    if(nrow(centers)==K) K <- NULL
  }
  if(!is.null(centers) && ncol(centers)!=ncol(x)) stop("If centers is provided, then ncol(centers) must equal ncol(x).")
  if(is.null(wbounds)) wbounds <- seq(1.1, sqrt(ncol(x)), len=20)
  if(min(wbounds)<=1) stop("wbounds should be greater than 1")
  wbounds <- c(wbounds) # In case wbounds is a single number, turn it into a vector
  out <- list()
  if(!is.null(K)) Cs <- kmeans(x, centers=K, nstart=nstart)$cluster
  if(is.null(K)) Cs <- kmeans(x, centers=centers)$cluster
  for(i in 1:length(wbounds)){
    if(length(wbounds)>1 && !silent) cat(i,fill=FALSE)
    ws <- rep(1/sqrt(ncol(x)), ncol(x)) # Start with equal weights on each feature
    ws.old <- rnorm(ncol(x))
    store.bcss.ws <- NULL
    niter <- 0
    while((sum(abs(ws-ws.old))/sum(abs(ws.old)))>1e-4 && niter<maxiter){
      if(!silent) cat(niter, fill=FALSE)
      niter <- niter+1
      ws.old <- ws
      if(!is.null(K)){
        if(niter>1) Cs <- UpdateCs(x, K, ws, Cs) # if niter=1, no need to update!!
      } else {
        if(niter>1) Cs <- UpdateCs(x, nrow(centers), ws, Cs) # if niter=1, no need to update!!
      }
      ws <- UpdateWs(x, Cs, wbounds[i])
      store.bcss.ws <- c(store.bcss.ws, sum(GetWCSS(x, Cs)$bcss.perfeature*ws))
    }
    out[[i]] <- list(ws=ws, Cs=Cs, wcss=GetWCSS(x, Cs, ws), crit=store.bcss.ws, wbound=wbounds[i])
  }
  if(!silent) cat(fill=TRUE)
  #  if(length(wbounds)==1){
  #    out <- out[[1]]
  #    class(out) <- "kmeanssparse"
  #    return(out)
  #  }
  #  class(out) <- "multikmeanssparse"
  #  return(out)
  class(out) <- "KmeansSparseCluster"
  return(out)
}

plot.KmeansSparseCluster <- function(x,...){
  if(length(x)>1){
    N <- length(x)
    par(mfrow=c(ceiling(N/2),2))
    for(i in 1:N){
      plot(x[[i]]$ws, main=paste("Wbound is ", sep="", round(x[[i]]$wbound,3)), xlab="Feature Index", ylab="Wj")
    }
  } else {
    x <- x[[1]]
    plot(x$ws, main=paste("Wbound is ", sep="", round(x$wbound,3)), xlab="Feature Index", ylab="Wj")
  }
}

#plot.kmeanssparse <- function(x,...){
#
#}

PrintIt <- function(x){
  cat("Number of non-zero weights: ", sum(x$ws!=0), fill=TRUE)
  cat("Sum of weights: ", sum(x$ws), fill=TRUE)
  cat("Clustering: ", x$Cs, fill=TRUE)
  cat(fill=TRUE)
}

print.kmeanssparse <- function(x,...){
  cat("Wbound is ", x$wbound, ":", fill=TRUE)
  PrintIt(x)
}

print.KmeansSparseCluster <- function(x,...){
  for(i in 1:length(x)){
    cat("Wbound is ", x[[i]]$wbound, ":", fill=TRUE)
    PrintIt(x[[i]])
  }
}


# this file contains the hidden functions for the sparcl package

GetWCSS <- function(x, Cs, ws=NULL){
  wcss.perfeature <- numeric(ncol(x))
  for(k in unique(Cs)){
    whichers <- (Cs==k)
    if(sum(whichers)>1) wcss.perfeature <- wcss.perfeature + apply(scale(x[whichers,],center=TRUE, scale=FALSE)^2, 2, sum)
  }
  bcss.perfeature <- apply(scale(x, center=TRUE, scale=FALSE)^2, 2, sum)-wcss.perfeature
  if(!is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), wcss.ws=sum(wcss.perfeature*ws),
                               bcss.perfeature=bcss.perfeature))
  if(is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), bcss.perfeature=bcss.perfeature))
}

UpdateCs <- function(x, K, ws, Cs){
  x <- x[,ws!=0]
  z <- sweep(x, 2, sqrt(ws[ws!=0]), "*")
  nrowz <- nrow(z)
  mus <- NULL
  if(!is.null(Cs)){
    for(k in unique(Cs)){
      if(sum(Cs==k)>1) mus <- rbind(mus, apply(z[Cs==k,],2,mean))
      if(sum(Cs==k)==1) mus <- rbind(mus, z[Cs==k,])
    }
  }
  if(is.null(mus)){
    km <- kmeans(z, centers=K, nstart=10)
  } else {
    distmat <- as.matrix(dist(rbind(z, mus)))[1:nrowz, (nrowz+1):(nrowz+K)]
    nearest <- apply(distmat, 1, which.min)
    if(length(unique(nearest))==K){
      km <- kmeans(z, centers=mus)
    } else {
      km <- kmeans(z, centers=K, nstart=10)
    }
  }
  return(km$cluster)
}

#distmat <- function(x){
#  return(dist(x))
#}

BinarySearch <- function(argu,sumabs){
  if(l2n(argu)==0 || sum(abs(argu/l2n(argu)))<=sumabs) return(0)
  lam1 <- 0
  lam2 <- max(abs(argu))-1e-5
  iter <- 1
  while(iter<=15 && (lam2-lam1)>(1e-4)){
    su <- soft(argu,(lam1+lam2)/2)
    if(sum(abs(su/l2n(su)))<sumabs){
      lam2 <- (lam1+lam2)/2
    } else {
      lam1 <- (lam1+lam2)/2
    }
    iter <- iter+1
  }
  return((lam1+lam2)/2)
}

soft <- function(x,d){
  return(sign(x)*pmax(0, abs(x)-d))
}

l2n <- function(vec){
  return(sqrt(sum(vec^2)))
}

GetUW <- function(ds, wbound,niter,uorth,silent){
  # uorth would be a $n^2 x k$ matrix containing $k$ previous
  # dissimilarity matrices found, if we want to do sparse comp clust
  # after having already done $k$ of these things
  # Example:
  # out <- HierarchicalSparseCluster(x,wbound=5)
  # out2 <- HierarchicalSparseCluster(x,wbound=5, uorth=out$u)
  # Then out2 contains a sparse complementary clustering
  p <- ncol(ds)
  w <- rep(1/p, p)*wbound
  iter <- 1
  if(!is.null(uorth)){
    if(sum(abs(uorth-t(uorth)))>1e-10) stop("Uorth must be symmetric!!!!!!!!!!")
    uorth <- matrix(uorth[lower.tri(uorth)],ncol=1)
    uorth <- uorth/sqrt(sum(uorth^2))
  }
  u <- rnorm(nrow(ds))
  w <- rep(1, ncol(ds))
  w.old <- rnorm(ncol(ds))
  #  u.old <- rnorm(nrow(ds))
  while(iter<=niter && (sum(abs(w.old-w))/sum(abs(w.old)))>1e-4){#(sum(abs(u.old-u))/sum(abs(u.old)))>1e-4){ # was 1e-3 and involved ws until 11.12.09
    if(!silent) cat(iter,fill=FALSE)
    #    u.old <- u
    if(iter>1) u <- ds[,argw>=lam]%*%matrix(w[argw>=lam],ncol=1)
    if(iter==1) u <- ds%*%matrix(w,ncol=1)
    if(!is.null(uorth)) u <- u - uorth%*%(t(uorth)%*%u)
    iter <- iter+1
    u <- u/l2n(u)
    w.old <- w
    argw <- matrix(pmax(matrix(u,nrow=1)%*%ds,0),ncol=1)
    lam <- BinarySearch(argw,wbound)
    w <- soft(argw,lam) # Don't need to normalize b/c this will be done soon
    w <- w/l2n(w)
  }
  u <-  ds[,argw>=lam]%*%matrix(w[argw>=lam],ncol=1)/sum(w)
  if(!is.null(uorth)) u <- u - uorth%*%(t(uorth)%*%u)
  u <- u/l2n(u)
  w <- w/l2n(w)
  crit <- sum(u*(ds%*%matrix(w,ncol=1)))
  u2 <- matrix(0,nrow=ceiling(sqrt(2*length(u))),ncol=ceiling(sqrt(2*length(u))))
  u2[lower.tri(u2)] <- u
  u <- as.matrix(as.dist(u2))/sqrt(2)
  return(list(u=u, w=w, crit=crit))
}


UpdateWs <- function(x, Cs, l1bound){
  wcss.perfeature <- GetWCSS(x, Cs)$wcss.perfeature
  tss.perfeature <- GetWCSS(x, rep(1, nrow(x)))$wcss.perfeature
  if(sum(tss.perfeature<wcss.perfeature)!=0){
    index_flation<-which(tss.perfeature<wcss.perfeature)
    tss.perfeature[index_flation]<-wcss.perfeature[index_flation]
  }
  lam <- BinarySearch(-wcss.perfeature+tss.perfeature, l1bound)
  ws.unscaled <- soft(-wcss.perfeature+tss.perfeature,lam)
  return(ws.unscaled/l2n(ws.unscaled))
}


output.cluster.files.fun <- function(x,out,outputfile.prefix,genenames=NULL,genedesc=NULL){
  p=ncol(x)
  n=nrow(x)
  geneids=dimnames(x)[[2]]
  samplenames=dimnames(x)[[1]]
  if(is.null(geneids)) geneids <- paste("Gene", 1:ncol(x))
  if(is.null(samplenames)) samplenames <- paste("Sample",1:nrow(x))
  if(is.null(genenames)){genenames=geneids}
  if(is.null(genedesc)){genedesc <- rep("", ncol(x))}

  xx=x[,out$ws!=0]
  geneids.subset=geneids[out$ws!=0]
  genenames.subset=genenames[out$ws!=0]
  genedesc.subset=genedesc[out$ws!=0]
  pp=ncol(xx)
  sample.order=out$hc$order
  samplenames.o=samplenames[sample.order]
  arrynames=paste("ARRY",as.character(sample.order),"X",sep="")
  feature.order=1:pp
  if(!is.null(out$hc.features)){feature.order=out$hc.features$order}
  xx.o=xx[sample.order,feature.order]
  geneids.subset.o=geneids.subset[feature.order]
  genenames.subset.o=genenames.subset[feature.order]
  genedesc.subset.o=genedesc.subset[feature.order]
  genex=paste("GENE",as.character(1:pp),"X",sep="")
  genex.o=genex[feature.order]
  arrynames.o=arrynames[sample.order]

  # output cdt
  file=paste(outputfile.prefix,".cdt",sep="")
  xprefix <- rbind(c("GID","UID","NAME","GWEIGHT",samplenames.o),c("AID","","","",arrynames))
  xbig <- matrix(NA, ncol=(n+4),nrow=pp)
  for(i in 1:pp){
    xbig[i,] <- c(genex.o[i],geneids.subset.o[i],genedesc.subset.o[i],"1",xx.o[,i]) # was xx
  }
  xbig <- rbind(xprefix,xbig)
  write.table(file=file,xbig,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

  #output atr file
  atr=out$hc$merge
  atr.new=atr
  for(i in 1:nrow(atr)){
    for(j in 1:2){
      if(atr[i,j]<0){atr.new[i,j]=paste("ARRY",as.character(abs(atr[i,j])),"X",sep="")}
      if(atr[i,j]>0){atr.new[i,j]=paste("NODE",as.character(abs(atr[i,j])),"X",sep="")}
    }}
  col1=paste("NODE",as.character(1:nrow(atr.new)),"X",sep="")
  atr.new=cbind(col1,atr.new,1-out$hc$height/2)
  output.matrix(atr.new, paste(outputfile.prefix,".atr",sep=""))

  if(!is.null(out$hc.features)){
    #output gtr file
    gtr=out$hc.features$merge
    gtr.new=gtr
    for(i in 1:nrow(gtr)){
      for(j in 1:2){
        if(gtr[i,j]<0){gtr.new[i,j]=paste("GENE",as.character(abs(gtr[i,j])),"X",sep="")}
        if(gtr[i,j]>0){gtr.new[i,j]=paste("NODE",as.character(abs(gtr[i,j])),"X",sep="")}
      }}
    col1=paste("NODE",as.character(1:nrow(gtr.new)),"X",sep="")
    gtr.new=cbind(col1,gtr.new,1-out$hc.features$height/2)
    output.matrix(gtr.new, paste(outputfile.prefix,".gtr",sep=""))
  }

  return()
}


output.matrix <- function(x,file){
  write.table(file=file,x,quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}



read.gct <- function(file) {
  if (is.character(file))
    if (file == "")
      file <- stdin()
    else {
      file <- file(file, "r")
      on.exit(close(file))
    }
    if (!inherits(file, "connection"))
      stop("argument `file' must be a character string or connection")

    # line 1 version
    version <- readLines(file, n=1)

    # line 2 dimensions
    dimensions <- scan(file, what=list("integer", "integer"), nmax=1, quiet=TRUE)
    rows <- dimensions[[1]]
    columns <- dimensions[[2]]
    # line 3 Name\tDescription\tSample names...
    column.names <- read.table(file, header=FALSE, quote='', nrows=1, sep="\t", fill=FALSE, comment.char='')
    column.names <-column.names[3:length(column.names)]


    if(length(column.names)!=columns) {
      stop(paste("Number of sample names", length(column.names), "not equal to the number of columns", columns, "."))
    }

    colClasses <- c(rep(c("character"), 2), rep(c("double"), columns))

    x <- read.table(file, header=FALSE, quote="", row.names=NULL, comment.char="", sep="\t", colClasses=colClasses, fill=FALSE)
    row.descriptions <- as.character(x[,2])
    data <- as.matrix(x[seq(from=3, to=dim(x)[2], by=1)])

    column.names <- column.names[!is.na(column.names)]

    colnames(data) <- column.names
    row.names(data) <- x[,1]
    return(list(row.descriptions=row.descriptions, data=data))
}


extract.prefix <- function(file){
  #  i=0
  #  while(substring(file,i,i)!="." & (i <nchar(file))) {i=i+1}
  #  if(i==nchar(file)){stop('Error in file name')}
  #  pre=substring(file,1,i-1)
  #  return(pre)
  tmp <- strsplit(file,"\\.")[[1]]
  return(paste(tmp[-length(tmp)],collapse="."))

}


LD<-function(x, K.try = 2:10,n.resample){
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
  if(max(result_data)<0.8){
    maxk<-1
  }else{
    maxk<-K.try[which.max(result_data)]
  }

  return(maxk)
}




S4<-function(x, K.try = 2:10,n.resample,trim=0.05,cut.off=0.8){
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

      #per<-1-trim
      stat = rev(sapply(1:(dim(o.mtx)[1]), stat.calc1))
      stat<-sort(stat,decreasing = T)
      #plot(stat)
      #title(paste("plot of S4.Iter","_K=",K,sep=""))
      #barplot(stat,width = 0.832)
      #title(paste("k=",K,sep=""))
      #abline(v=n*0.95,col="red",lwd=3.5)
      U.min = mean(stat)
      return(U.min)
    })
  }
  result_data<-result_k(x,K.try)


  #if(specificity==TRUE){
  #  result<-apply(result_data,2,mean)
  #}else{
  #  result<-result_data[1,]
  #}
  #if((res.null==1)&(result_data==1)){
  #  maxk<-K.try[max(which(result_data == max(result_data)))]
  #}else{
  if(length(which(is.na(result_data)))!=0){
    result_data[which(is.na(result_data))]<-0
  }
  if(max(result_data)<cut.off){
    maxk<-1
  }else{
    maxk<-K.try[max(which(result_data == max(result_data)))]
  }

  #}
  res<-list(maxk=maxk,score=result_data)
  return(res)
}


cluster_CH<-function(x,K.try=2:10){
  km <- sapply(K.try,function(k){
    clus<-kmeans(x,k,nstart=100)
    return(calinhara(x,clus$cluster))})
  index<-which.max(km)
  return(K.try[index])
}

jump<-function(data=NULL,K=10,y=NULL,plotjumps=F,rand=100,fits=NULL,B=0,dist=NULL,trace=F){
  compute.jump <- function(dist,y,plotjumps=F,printresults=T){
    K <- length(dist)
    numb.y <- length(y)
    numbclust <- rep(0,numb.y)
    transdist <- matrix(0,numb.y,K+1)
    jumps <- matrix(0,numb.y,K)
    if (plotjumps)
      par(mfrow=c(numb.y,3))
    for (i in 1:numb.y){
      # Compute the transformed distortion
      transdist[i,] <- c(0,dist^(-y[i]))
      # Compute the jumps in transformed distortion
      jumps[i,] <- diff(transdist[i,])
      # Compute the maximum jump
      numbclust[i] <- order(-jumps[i,])[1]
      # Plot distortion, transformed distortion and jumps
      if (plotjumps){
        plot(1:K,dist,type='l',
             xlab="Number of Clusters",ylab="Distortion",main=paste("Y = ",y[i]))
        plot(0:K,transdist[i,],type='l',
             xlab="Number of Clusters",ylab="Transformed Distortion",main=paste("Y = ",y[i]))
        plot(1:K,jumps[i,],type='l',
             xlab="Number of Clusters",ylab="Jumps",main=paste("Y = ",y[i]))
        # Plot line and point to indicate maximum jump
        lines(rep(numbclust[i],2),c(0,jumps[i,numbclust[i]]),lty=3,lwd=3)
        points(numbclust[i],jumps[i,numbclust[i]],col=2,pch=19,cex=1.5)}
      # Report maximum jump
      if (printresults)
        print(paste("The maximum jump occurred at ",numbclust[i], "clusters with Y=",y[i]))
    }
    list(maxjump=numbclust,dist=dist,transdist=transdist[,-1],jumps=jumps)}

  kmeans.rndstart <- function(x, K, rand = 10)
  {
    fits <- sum((t(x) - apply(x, 2, mean))^2)
    iter.max <- 10
    # Run kmeans for 2 to K clusters
    for (k in 2:K){
      Z=kmeans(x,k,nstart=rand)
      fits=c(fits,sum(Z$withinss))
    }
    fits
  }



  if (!is.null(data)){
    # Compute the kmeans fit to the data
    if (is.null(fits))
      fits <- kmeans.rndstart(data,K,rand)
    if (is.null(y))
      y <- dim(data)[2]/2
    n <- nrow(data)
    p <- ncol(data)
    # Compute the distortion associated with the kmeans fit
    dist<- fits/(n*p)
  }
  # Call the compute.jump function to produce plots and calculate
  # maximum jump for each value of Y
  jump.results <- compute.jump(dist,y,plotjumps)
  jump.results$fits <- fits
  # Implement bootstrap routine
  if (B>0 & !is.null(data)){
    n <- nrow(data)
    boot.results <- matrix(0,length(y),K)
    bootdist <- rep(0,K)
    for (b in 1:B){
      if (trace)
        print(paste("Bootstrap Iteration ",b))
      # Make bootstrap data
      bootdata <- data[sample(1:n,replace=T),]
      # Get kmeans fit to the bootstrap data
      bootfits <- kmeans.rndstart(bootdata,K,rand)
      # Compute bootstrap distortion and maximum jumps
      for (k in 1:K)
        bootdist[k] <- sum(bootfits[[k]]$within)/(n*p)
      bootmaxjumps <- compute.jump(bootdist,y,plotjumps=F,printresults=F)$maxjump
      for (j in 1:length(y))
        boot.results[j,bootmaxjumps[j]] <-  boot.results[j,bootmaxjumps[j]]+1
    }
    # Calculate proportions of each number of clusters chosen
    jump.results$boot.result <- round(boot.results/B,3)
    #    for (j in 1:length(y))
    #      print(paste(jump.results$boot.result[j,jump.results$maxjump[j]]*100,"% of bootstrap iterations corresponding to ",jump.results$maxjump[j], "clusters with Y=",y[j]))
  }
  jump.results}

prediction.strength <- function (xdata, Gmin = 2, Gmax = 10, M = 50, clustermethod = kmeans,
                                 classification = "centroid", cutoff = 0.8, nnk = 1, distances = inherits(xdata,
                                                                                                          "dist"), count = FALSE, ...)
{
  xdata <- as.matrix(xdata)
  n <- nrow(xdata)
  nf <- c(floor(n/2), n - floor(n/2))
  indvec <- clcenters <- clusterings <- jclusterings <- classifications <- list()
  corrpred <- list()
  for (k in Gmin:Gmax) {
    if (count)
      cat(k, " clusters\n")
    corrpred[[k]] <- numeric(0)
    for (l in 1:M) {
      nperm <- sample(n, n)
      if (count)
        cat(" Run ", l, "\n")
      indvec[[l]] <- list()
      indvec[[l]][[1]] <- nperm[1:nf[1]]
      indvec[[l]][[2]] <- nperm[(nf[1] + 1):n]
      for (i in 1:2) {
        if (distances)
          clusterings[[i]] <- kmeans(as.dist(xdata[indvec[[l]][[i]],]),centers = k, nstart = 100,...)
        else clusterings[[i]] <- kmeans(xdata[indvec[[l]][[i]],],centers = k, nstart = 100,...)

        jclusterings[[i]] <- rep(-1, n)
        jclusterings[[i]][indvec[[l]][[i]]] <- clusterings[[i]]$cluster
        centroids <- NULL
        if (classification == "centroid") {
          if (identical(clustermethod, kmeans))
            centroids <- clusterings[[i]]$centers
          if (identical(clustermethod, claraCBI))
            centroids <- clusterings[[i]]$result$medoids
        }
        j <- 3 - i
        if (distances)
          classifications[[j]] <- classifdist(as.dist(xdata),
                                              jclusterings[[i]], method = classification,
                                              centroids = centroids, nnk = nnk)[indvec[[l]][[j]]]
        else classifications[[j]] <- classifnp(xdata,
                                               jclusterings[[i]], method = classification,
                                               centroids = centroids, nnk = nnk)[indvec[[l]][[j]]]
      }
      ps <- matrix(0, nrow = 2, ncol = k)
      for (i in 1:2) {
        ctable <- table(clusterings[[i]]$cluster, classifications[[i]])
        for (kk in 1:k) {
          ps[i, kk] <- sum(ctable[kk, ]^2 - ctable[kk,
                                                   ])
          cpik <- clusterings[[i]]$cluster == kk
          nik <- sum(cpik)
          if (nik > 1)
            ps[i, kk] <- ps[i, kk]/(nik * (nik - 1))
          else ps[i, kk] <- 1
        }
      }
      corrpred[[k]][l] <- mean(c(min(ps[1, ]), min(ps[2,
                                                      ])))
    }
  }
  mean.pred <- numeric(0)
  if (Gmin > 1)
    mean.pred <- c(1)
  if (Gmin > 2)
    mean.pred <- c(mean.pred, rep(NA, Gmin - 2))
  for (k in Gmin:Gmax) mean.pred <- c(mean.pred, mean(corrpred[[k]]))
  optimalk <- max(which(mean.pred > cutoff))
  out <- list(predcorr = corrpred, mean.pred = mean.pred, optimalk = optimalk,
              cutoff = cutoff, method = clusterings[[1]]$clustermethod,
              Gmax = Gmax, M = M)
  class(out) <- "predstr"
  out
}








hartigan1<-function(x,min.nc = 2, max.nc = 10){
  n<-nrow(x)
  p<-ncol(x)
  # between within sum of squares
  bss<-0
  wss <- (nrow(x)-1)*sum(apply(x,2,var))
  # four statistics when k=1

  # between and within ss variying k from 2 to 11
  for (k in 2:11) {
    kmc<-kmeans(x,centers=k,nstart=20)
    wss[k] <-kmc$tot.withinss
    bss[k] <-kmc$betweenss
    # silhouette
  }

  H<-(wss[1]/wss[2]-1)*(n-2)

  # four statistics varying k from 2 to 10
  for(k in 2:10){

    # Hartigan
    H[k]<-(wss[k]/wss[k+1] -1)*(n-k-1)

  }

  if(sum(which(H<=10))==0){
    H.k<-10
  } else{
    H.k<-min(which(H<=10))
  }
  return(H.k)
}


NBClust1<-function (data = NULL, diss = NULL, distance = "euclidean",
                    min.nc = 2, max.nc = 15, method = NULL, index = "all", alphaBeale = 0.1)
{
  x <- 0
  min_nc <- min.nc
  max_nc <- max.nc
  if (is.null(method))
    stop("method is NULL")
  method <- pmatch(method, c("ward.D2", "single", "complete",
                             "average", "mcquitty", "median", "centroid", "kmeans",
                             "ward.D"))
  indice <- pmatch(index, c("kl", "ch", "hartigan", "ccc",
                            "scott", "marriot", "trcovw", "tracew", "friedman",
                            "rubin", "cindex", "db", "silhouette", "duda", "pseudot2",
                            "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey",
                            "mcclain", "gamma", "gplus", "tau", "dunn", "hubert",
                            "sdindex", "dindex", "sdbw", "all", "alllong"))
  if (is.na(indice))
    stop("invalid clustering index")
  if (indice == -1)
    stop("ambiguous index")
  if ((indice == 3) || (indice == 5) || (indice == 6) || (indice ==
                                                          7) || (indice == 8) || (indice == 9) || (indice == 10) ||
      (indice == 11) || (indice == 18) || (indice == 27) ||
      (indice == 29) || (indice == 31) || (indice == 32)) {
    if ((max.nc - min.nc) < 2)
      stop("The difference between the minimum and the maximum number of clusters must be at least equal to 2")
  }
  if (is.null(data)) {
    if (method == 8) {
      stop("\n", "method = kmeans, data matrix is needed")
    }
    else {
      if ((indice == 1) || (indice == 2) || (indice ==
                                             3) || (indice == 4) || (indice == 5) || (indice ==
                                                                                      6) || (indice == 7) || (indice == 8) || (indice ==
                                                                                                                               9) || (indice == 10) || (indice == 12) || (indice ==
                                                                                                                                                                          14) || (indice == 15) || (indice == 16) || (indice ==
                                                                                                                                                                                                                      17) || (indice == 18) || (indice == 19) || (indice ==
                                                                                                                                                                                                                                                                  20) || (indice == 23) || (indice == 24) || (indice ==
                                                                                                                                                                                                                                                                                                              25) || (indice == 27) || (indice == 28) || (indice ==
                                                                                                                                                                                                                                                                                                                                                          29) || (indice == 30) || (indice == 31) || (indice ==
                                                                                                                                                                                                                                                                                                                                                                                                      32))
        stop("\n", "Data matrix is needed. Only frey, mcclain, cindex, sihouette and dunn can be computed.",
             "\n")
      if (is.null(diss))
        stop("data matrix and dissimilarity matrix are both null")
      else cat("\n", "Only frey, mcclain, cindex, sihouette and dunn can be computed. To compute the other indices, data matrix is needed",
               "\n")
    }
  }
  else {
    jeu1 <- as.matrix(data)
    numberObsBefore <- dim(jeu1)[1]
    jeu <- na.omit(jeu1)
    nn <- numberObsAfter <- dim(jeu)[1]
    pp <- dim(jeu)[2]
    TT <- t(jeu) %*% jeu
    sizeEigenTT <- length(eigen(TT)$value)
    eigenValues <- eigen(TT/(nn - 1))$value
    if (any(indice == 4) || (indice == 5) || (indice ==
                                              6) || (indice == 7) || (indice == 8) || (indice ==
                                                                                       9) || (indice == 10) || (indice == 31) || (indice ==
                                                                                                                                  32)) {
      for (i in 1:sizeEigenTT) {
        if (eigenValues[i] < 0) {
          stop("The TSS matrix is indefinite. There must be too many missing values. The index cannot be calculated.")
        }
      }
      s1 <- sqrt(eigenValues)
      ss <- rep(1, sizeEigenTT)
      for (i in 1:sizeEigenTT) {
        if (s1[i] != 0)
          ss[i] = s1[i]
      }
      vv <- prod(ss)
    }
  }
  if (is.null(distance))
    distanceM <- 7
  if (!is.null(distance))
    distanceM <- pmatch(distance, c("euclidean", "maximum",
                                    "manhattan", "canberra", "binary", "minkowski"))
  if (is.na(distanceM)) {
    stop("invalid distance")
  }
  if (is.null(diss)) {
    if (distanceM == 1) {
      md <- dist(jeu, method = "euclidean")
    }
    if (distanceM == 2) {
      md <- dist(jeu, method = "maximum")
    }
    if (distanceM == 3) {
      md <- dist(jeu, method = "manhattan")
    }
    if (distanceM == 4) {
      md <- dist(jeu, method = "canberra")
    }
    if (distanceM == 5) {
      md <- dist(jeu, method = "binary")
    }
    if (distanceM == 6) {
      md <- dist(jeu, method = "minkowski")
    }
    if (distanceM == 7) {
      stop("dissimilarity matrix and distance are both NULL")
    }
  }
  if (!is.null(diss)) {
    if ((distanceM == 1) || (distanceM == 2) || (distanceM ==
                                                 3) || (distanceM == 4) || (distanceM == 5) || (distanceM ==
                                                                                                6))
      stop("dissimilarity matrix and distance are both not null")
    else md <- diss
  }
  res <- array(0, c(max_nc - min_nc + 1, 30))
  x_axis <- min_nc:max_nc
  resCritical <- array(0, c(max_nc - min_nc + 1, 4))
  rownames(resCritical) <- min_nc:max_nc
  colnames(resCritical) <- c("CritValue_Duda", "CritValue_PseudoT2",
                             "Fvalue_Beale", "CritValue_Gap")
  rownames(res) <- min_nc:max_nc
  colnames(res) <- c("KL", "CH", "Hartigan", "CCC", "Scott",
                     "Marriot", "TrCovW", "TraceW", "Friedman", "Rubin",
                     "Cindex", "DB", "Silhouette", "Duda", "Pseudot2", "Beale",
                     "Ratkowsky", "Ball", "Ptbiserial", "Gap", "Frey", "McClain",
                     "Gamma", "Gplus", "Tau", "Dunn", "Hubert", "SDindex",
                     "Dindex", "SDbw")
  if (is.na(method))
    stop("invalid clustering method")
  if (method == -1)
    stop("ambiguous method")
  if (method == 1) {
    hc <- hclust(md, method = "ward.D2")
  }
  if (method == 2) {
    hc <- hclust(md, method = "single")
  }
  if (method == 3) {
    hc <- hclust(md, method = "complete")
  }
  if (method == 4) {
    hc <- hclust(md, method = "average")
  }
  if (method == 5) {
    hc <- hclust(md, method = "mcquitty")
  }
  if (method == 6) {
    hc <- hclust(md, method = "median")
  }
  if (method == 7) {
    hc <- hclust(md, method = "centroid")
  }
  if (method == 9) {
    hc <- hclust(md, method = "ward.D")
  }
  centers <- function(cl, x) {
    x <- as.matrix(x)
    n <- length(cl)
    k <- max(cl)
    centers <- matrix(nrow = k, ncol = ncol(x))
    {
      for (i in 1:k) {
        for (j in 1:ncol(x)) {
          centers[i, j] <- mean(x[cl == i, j])
        }
      }
    }
    return(centers)
  }
  Average.scattering <- function(cl, x) {
    x <- as.matrix(x)
    n <- length(cl)
    k <- max(cl)
    centers.matrix <- centers(cl, x)
    cluster.size <- numeric(0)
    variance.clusters <- matrix(0, ncol = ncol(x), nrow = k)
    var <- matrix(0, ncol = ncol(x), nrow = k)
    for (u in 1:k) cluster.size[u] <- sum(cl == u)
    for (u in 1:k) {
      for (j in 1:ncol(x)) {
        for (i in 1:n) {
          if (cl[i] == u)
            variance.clusters[u, j] <- variance.clusters[u,
                                                         j] + (x[i, j] - centers.matrix[u, j])^2
        }
      }
    }
    for (u in 1:k) {
      for (j in 1:ncol(x)) variance.clusters[u, j] = variance.clusters[u,
                                                                       j]/cluster.size[u]
    }
    variance.matrix <- numeric(0)
    for (j in 1:ncol(x)) variance.matrix[j] = var(x[, j]) *
      (n - 1)/n
    Somme.variance.clusters <- 0
    for (u in 1:k) Somme.variance.clusters <- Somme.variance.clusters +
      sqrt((variance.clusters[u, ] %*% (variance.clusters[u,
                                                          ])))
    stdev <- (1/k) * sqrt(Somme.variance.clusters)
    scat <- (1/k) * (Somme.variance.clusters/sqrt(variance.matrix %*%
                                                    variance.matrix))
    scat <- list(stdev = stdev, centers = centers.matrix,
                 variance.intraclusters = variance.clusters, scatt = scat)
    return(scat)
  }
  density.clusters <- function(cl, x) {
    x <- as.matrix(x)
    k <- max(cl)
    n <- length(cl)
    distance <- matrix(0, ncol = 1, nrow = n)
    density <- matrix(0, ncol = 1, nrow = k)
    centers.matrix <- centers(cl, x)
    stdev <- Average.scattering(cl, x)$stdev
    for (i in 1:n) {
      u = 1
      while (cl[i] != u) u <- u + 1
      for (j in 1:ncol(x)) {
        distance[i] <- distance[i] + (x[i, j] - centers.matrix[u,
                                                               j])^2
      }
      distance[i] <- sqrt(distance[i])
      if (distance[i] <= stdev)
        density[u] = density[u] + 1
    }
    dens <- list(distance = distance, density = density)
    return(dens)
  }
  density.bw <- function(cl, x) {
    x <- as.matrix(x)
    k <- max(cl)
    n <- length(cl)
    centers.matrix <- centers(cl, x)
    stdev <- Average.scattering(cl, x)$stdev
    density.bw <- matrix(0, ncol = k, nrow = k)
    u <- 1
    for (u in 1:k) {
      for (v in 1:k) {
        if (v != u) {
          distance <- matrix(0, ncol = 1, nrow = n)
          moy <- (centers.matrix[u, ] + centers.matrix[v,
                                                       ])/2
          for (i in 1:n) {
            if ((cl[i] == u) || (cl[i] == v)) {
              for (j in 1:ncol(x)) {
                distance[i] <- distance[i] + (x[i, j] -
                                                moy[j])^2
              }
              distance[i] <- sqrt(distance[i])
              if (distance[i] <= stdev) {
                density.bw[u, v] <- density.bw[u, v] +
                  1
              }
            }
          }
        }
      }
    }
    density.clust <- density.clusters(cl, x)$density
    S <- 0
    for (u in 1:k) for (v in 1:k) {
      if (max(density.clust[u], density.clust[v]) != 0)
        S = S + (density.bw[u, v]/max(density.clust[u],
                                      density.clust[v]))
    }
    density.bw <- S/(k * (k - 1))
    return(density.bw)
  }
  Dis <- function(cl, x) {
    x <- as.matrix(x)
    k <- max(cl)
    centers.matrix <- centers(cl, x)
    Distance.centers <- dist(centers.matrix)
    Dmin <- min(Distance.centers)
    Dmax <- max(Distance.centers)
    Distance.centers <- as.matrix(Distance.centers)
    s2 <- 0
    for (u in 1:k) {
      s1 = 0
      for (j in 1:ncol(Distance.centers)) {
        s1 <- s1 + Distance.centers[u, j]
      }
      s2 <- s2 + 1/s1
    }
    Dis <- (Dmax/Dmin) * s2
    return(Dis)
  }
  Index.Hubert <- function(x, cl) {
    k <- max(cl)
    n <- dim(x)[1]
    y <- matrix(0, ncol = dim(x)[2], nrow = n)
    P <- as.matrix(md)
    meanP <- mean(P)
    variance.matrix <- numeric(0)
    md <- dist(x, method = "euclidean")
    for (j in 1:n) variance.matrix[j] = var(P[, j]) * (n -
                                                         1)/n
    varP <- sqrt(variance.matrix %*% variance.matrix)
    centers.clusters <- centers(cl, x)
    for (i in 1:n) {
      for (u in 1:k) {
        if (cl[i] == u)
          y[i, ] <- centers.clusters[u, ]
      }
    }
    Q <- as.matrix(dist(y, method = "euclidean"))
    meanQ <- mean(Q)
    for (j in 1:n) variance.matrix[j] = var(Q[, j]) * (n -
                                                         1)/n
    varQ <- sqrt(variance.matrix %*% variance.matrix)
    M <- n * (n - 1)/2
    S <- 0
    n1 <- n - 1
    for (i in 1:n1) {
      j <- i + 1
      while (j <= n) {
        S <- S + (P[i, j] - meanP) * (Q[i, j] - meanQ)
        j <- j + 1
      }
    }
    gamma <- S/(M * varP * varQ)
    return(gamma)
  }
  Index.sPlussMoins <- function(cl1, md) {
    cn1 <- max(cl1)
    n1 <- length(cl1)
    dmat <- as.matrix(md)
    average.distance <- median.distance <- separation <- cluster.size <- within.dist1 <- between.dist1 <- numeric(0)
    separation.matrix <- matrix(0, ncol = cn1, nrow = cn1)
    di <- list()
    for (u in 1:cn1) {
      cluster.size[u] <- sum(cl1 == u)
      du <- as.dist(dmat[cl1 == u, cl1 == u])
      within.dist1 <- c(within.dist1, du)
      average.distance[u] <- mean(du)
      median.distance[u] <- median(du)
      bv <- numeric(0)
      for (v in 1:cn1) {
        if (v != u) {
          suv <- dmat[cl1 == u, cl1 == v]
          bv <- c(bv, suv)
          if (u < v) {
            separation.matrix[u, v] <- separation.matrix[v,
                                                         u] <- min(suv)
            between.dist1 <- c(between.dist1, suv)
          }
        }
      }
    }
    nwithin1 <- length(within.dist1)
    nbetween1 <- length(between.dist1)
    meanwithin1 <- mean(within.dist1)
    meanbetween1 <- mean(between.dist1)
    s.plus <- s.moins <- 0
    for (k in 1:nwithin1) {
      s.plus <- s.plus + (colSums(outer(between.dist1,
                                        within.dist1[k], ">")))
      s.moins <- s.moins + (colSums(outer(between.dist1,
                                          within.dist1[k], "<")))
    }
    Index.Gamma <- (s.plus - s.moins)/(s.plus + s.moins)
    Index.Gplus <- (2 * s.moins)/(n1 * (n1 - 1))
    t.tau <- (nwithin1 * nbetween1) - (s.plus + s.moins)
    Index.Tau <- (s.plus - s.moins)/(((n1 * (n1 - 1)/2 -
                                         t.tau) * (n1 * (n1 - 1)/2))^(1/2))
    results <- list(gamma = Index.Gamma, gplus = Index.Gplus,
                    tau = Index.Tau)
    return(results)
  }
  Index.15and28 <- function(cl1, cl2, md) {
    cn1 <- max(cl1)
    n1 <- length(cl1)
    dmat <- as.matrix(md)
    average.distance <- median.distance <- separation <- cluster.size <- within.dist1 <- between.dist1 <- numeric(0)
    separation.matrix <- matrix(0, ncol = cn1, nrow = cn1)
    di <- list()
    for (u in 1:cn1) {
      cluster.size[u] <- sum(cl1 == u)
      du <- as.dist(dmat[cl1 == u, cl1 == u])
      within.dist1 <- c(within.dist1, du)
      for (v in 1:cn1) {
        if (v != u) {
          suv <- dmat[cl1 == u, cl1 == v]
          if (u < v) {
            separation.matrix[u, v] <- separation.matrix[v,
                                                         u] <- min(suv)
            between.dist1 <- c(between.dist1, suv)
          }
        }
      }
    }
    cn2 <- max(cl2)
    n2 <- length(cl2)
    dmat <- as.matrix(md)
    average.distance <- median.distance <- separation <- cluster.size <- within.dist2 <- between.dist2 <- numeric(0)
    separation.matrix <- matrix(0, ncol = cn2, nrow = cn2)
    di <- list()
    for (w in 1:cn2) {
      cluster.size[w] <- sum(cl2 == w)
      dw <- as.dist(dmat[cl2 == w, cl2 == w])
      within.dist2 <- c(within.dist2, dw)
      bx <- numeric(0)
      for (x in 1:cn2) {
        if (x != w) {
          swx <- dmat[cl2 == w, cl2 == x]
          bx <- c(bx, swx)
          if (w < x) {
            separation.matrix[w, x] <- separation.matrix[x,
                                                         w] <- min(swx)
            between.dist2 <- c(between.dist2, swx)
          }
        }
      }
    }
    nwithin1 <- length(within.dist1)
    nbetween1 <- length(between.dist1)
    meanwithin1 <- mean(within.dist1)
    meanbetween1 <- mean(between.dist1)
    meanwithin2 <- mean(within.dist2)
    meanbetween2 <- mean(between.dist2)
    Index.15 <- (meanbetween2 - meanbetween1)/(meanwithin2 -
                                                 meanwithin1)
    Index.28 <- (meanwithin1/nwithin1)/(meanbetween1/nbetween1)
    results <- list(frey = Index.15, mcclain = Index.28)
    return(results)
  }
  Indice.ptbiserial <- function(x, md, cl1) {
    nn <- dim(x)[1]
    pp <- dim(x)[2]
    md2 <- as.matrix(md)
    m01 <- array(NA, c(nn, nn))
    nbr <- (nn * (nn - 1))/2
    pb <- array(0, c(nbr, 2))
    m3 <- 1
    for (m1 in 2:nn) {
      m12 <- m1 - 1
      for (m2 in 1:m12) {
        if (cl1[m1] == cl1[m2])
          m01[m1, m2] <- 0
        if (cl1[m1] != cl1[m2])
          m01[m1, m2] <- 1
        pb[m3, 1] <- m01[m1, m2]
        pb[m3, 2] <- md2[m1, m2]
        m3 <- m3 + 1
      }
    }
    y <- pb[, 1]
    x <- pb[, 2]
    biserial.cor <- function(x, y, use = c("all.obs", "complete.obs"),
                             level = 1) {
      if (!is.numeric(x))
        stop("'x' must be a numeric variable.\n")
      y <- as.factor(y)
      if (length(levs <- levels(y)) > 2)
        stop("'y' must be a dichotomous variable.\n")
      if (length(x) != length(y))
        stop("'x' and 'y' do not have the same length")
      use <- match.arg(use)
      if (use == "complete.obs") {
        cc.ind <- complete.cases(x, y)
        x <- x[cc.ind]
        y <- y[cc.ind]
      }
      ind <- y == levs[level]
      diff.mu <- mean(x[ind]) - mean(x[!ind])
      prob <- mean(ind)
      diff.mu * sqrt(prob * (1 - prob))/sd(x)
    }
    ptbiserial <- biserial.cor(x = pb[, 2], y = pb[, 1],
                               level = 2)
    return(ptbiserial)
  }
  Indices.WKWL <- function(x, cl1 = cl1, cl2 = cl2) {
    dim2 <- dim(x)[2]
    wss <- function(x) {
      x <- as.matrix(x)
      n <- length(x)
      centers <- matrix(nrow = 1, ncol = ncol(x))
      if (ncol(x) == 1) {
        centers[1, ] <- mean(x)
      }
      if (is.null(dim(x))) {
        bb <- matrix(x, byrow = FALSE, nrow = 1, ncol = ncol(x))
        centers[1, ] <- apply(bb, 2, mean)
      }
      else {
        centers[1, ] <- apply(x, 2, mean)
      }
      x.2 <- sweep(x, 2, centers[1, ], "-")
      withins <- sum(x.2^2)
      wss <- sum(withins)
      return(wss)
    }
    ncg1 <- 1
    ncg1max <- max(cl1)
    while ((sum(cl1 == ncg1) == sum(cl2 == ncg1)) && ncg1 <=
           ncg1max) {
      ncg1 <- ncg1 + 1
    }
    g1 <- ncg1
    ncg2 <- max(cl2)
    nc2g2 <- ncg2 - 1
    while ((sum(cl1 == nc2g2) == sum(cl2 == ncg2)) && nc2g2 >=
           1) {
      ncg2 <- ncg2 - 1
      nc2g2 <- nc2g2 - 1
    }
    g2 <- ncg2
    NK <- sum(cl2 == g1)
    WK.x <- x[cl2 == g1, ]
    WK <- wss(x = WK.x)
    NL <- sum(cl2 == g2)
    WL.x <- x[cl2 == g2, ]
    WL <- wss(x = WL.x)
    NM <- sum(cl1 == g1)
    WM.x <- x[cl1 == g1, ]
    WM <- wss(x = WM.x)
    duda <- (WK + WL)/WM
    BKL <- WM - WK - WL
    pseudot2 <- BKL/((WK + WL)/(NK + NL - 2))
    beale <- (BKL/(WK + WL))/(((NM - 1)/(NM - 2)) * (2^(2/dim2) -
                                                       1))
    results <- list(duda = duda, pseudot2 = pseudot2, NM = NM,
                    NK = NK, NL = NL, beale = beale)
    return(results)
  }
  Indices.WBT <- function(x, cl, P, s, vv) {
    n <- dim(x)[1]
    pp <- dim(x)[2]
    qq <- max(cl)
    z <- matrix(0, ncol = qq, nrow = n)
    clX <- as.matrix(cl)
    for (i in 1:n) for (j in 1:qq) {
      z[i, j] == 0
      if (clX[i, 1] == j) {
        z[i, j] = 1
      }
    }
    xbar <- solve(t(z) %*% z) %*% t(z) %*% x
    B <- t(xbar) %*% t(z) %*% z %*% xbar
    W <- P - B
    marriot <- (qq^2) * det(W)
    trcovw <- sum(diag(cov(W)))
    tracew <- sum(diag(W))
    if (det(W) != 0)
      scott <- n * log(det(P)/det(W))
    else {
      cat("Error: division by zero!")
    }
    friedman <- sum(diag(solve(W) * B))
    rubin <- sum(diag(P))/sum(diag(W))
    R2 <- 1 - sum(diag(W))/sum(diag(P))
    v1 <- 1
    u <- rep(0, pp)
    c <- (vv/(qq))^(1/pp)
    u <- s/c
    k1 <- sum((u >= 1) == TRUE)
    p1 <- min(k1, qq - 1)
    if (all(p1 > 0, p1 < pp)) {
      for (i in 1:p1) v1 <- v1 * s[i]
      c <- (v1/(qq))^(1/p1)
      u <- s/c
      b1 <- sum(1/(n + u[1:p1]))
      b2 <- sum(u[p1 + 1:pp]^2/(n + u[p1 + 1:pp]), na.rm = TRUE)
      E_R2 <- 1 - ((b1 + b2)/sum(u^2)) * ((n - qq)^2/n) *
        (1 + 4/n)
      ccc <- log((1 - E_R2)/(1 - R2)) * (sqrt(n * p1/2)/((0.001 +
                                                            E_R2)^1.2))
    }
    else {
      b1 <- sum(1/(n + u))
      E_R2 <- 1 - (b1/sum(u^2)) * ((n - qq)^2/n) * (1 +
                                                      4/n)
      ccc <- log((1 - E_R2)/(1 - R2)) * (sqrt(n * pp/2)/((0.001 +
                                                            E_R2)^1.2))
    }
    results <- list(ccc = ccc, scott = scott, marriot = marriot,
                    trcovw = trcovw, tracew = tracew, friedman = friedman,
                    rubin = rubin)
    return(results)
  }
  Indices.Traces <- function(data, d, clall, index = "all") {
    x <- data
    cl0 <- clall[, 1]
    cl1 <- clall[, 2]
    cl2 <- clall[, 3]
    clall <- clall
    nb.cl0 <- table(cl0)
    nb.cl1 <- table(cl1)
    nb.cl2 <- table(cl2)
    nb1.cl0 <- sum(nb.cl0 == 1)
    nb1.cl1 <- sum(nb.cl1 == 1)
    nb1.cl2 <- sum(nb.cl2 == 1)
    gss <- function(x, cl, d) {
      n <- length(cl)
      k <- max(cl)
      centers <- matrix(nrow = k, ncol = ncol(x))
      for (i in 1:k) {
        if (ncol(x) == 1) {
          centers[i, ] <- mean(x[cl == i, ])
        }
        if (is.null(dim(x[cl == i, ]))) {
          bb <- matrix(x[cl == i, ], byrow = FALSE,
                       nrow = 1, ncol = ncol(x))
          centers[i, ] <- apply(bb, 2, mean)
        }
        else {
          centers[i, ] <- apply(x[cl == i, ], 2, mean)
        }
      }
      allmean <- apply(x, 2, mean)
      dmean <- sweep(x, 2, allmean, "-")
      allmeandist <- sum(dmean^2)
      withins <- rep(0, k)
      x.2 <- (x - centers[cl, ])^2
      for (i in 1:k) {
        withins[i] <- sum(x.2[cl == i, ])
      }
      wgss <- sum(withins)
      bgss <- allmeandist - wgss
      results <- list(wgss = wgss, bgss = bgss, centers = centers)
      return(results)
    }
    vargss <- function(x, clsize, varwithins) {
      nvar <- dim(x)[2]
      n <- sum(clsize)
      k <- length(clsize)
      varallmean <- rep(0, nvar)
      varallmeandist <- rep(0, nvar)
      varwgss <- rep(0, nvar)
      for (l in 1:nvar) varallmean[l] <- mean(x[, l])
      vardmean <- sweep(x, 2, varallmean, "-")
      for (l in 1:nvar) {
        varallmeandist[l] <- sum((vardmean[, l])^2)
        varwgss[l] <- sum(varwithins[, l])
      }
      varbgss <- varallmeandist - varwgss
      vartss <- varbgss + varwgss
      zvargss <- list(vartss = vartss, varbgss = varbgss)
      return(zvargss)
    }
    varwithinss <- function(x, centers, cluster) {
      nrow <- dim(centers)[1]
      nvar <- dim(x)[2]
      varwithins <- matrix(0, nrow, nvar)
      x <- (x - centers[cluster, ])^2
      for (l in 1:nvar) {
        for (k in 1:nrow) {
          varwithins[k, l] <- sum(x[cluster == k, l])
        }
      }
      return(varwithins)
    }
    indice.kl <- function(x, clall, d = NULL, centrotypes = "centroids") {
      if (nb1.cl1 > 0) {
        KL <- NA
      }
      if (sum(c("centroids", "medoids") == centrotypes) ==
          0)
        stop("Wrong centrotypes argument")
      if ("medoids" == centrotypes && is.null(d))
        stop("For argument centrotypes = 'medoids' d cannot be null")
      if (!is.null(d)) {
        if (!is.matrix(d)) {
          d <- as.matrix(d)
        }
        row.names(d) <- row.names(x)
      }
      m <- ncol(x)
      g <- k <- max(clall[, 2])
      KL <- abs((g - 1)^(2/m) * gss(x, clall[, 1], d)$wgss -
                  g^(2/m) * gss(x, clall[, 2], d)$wgss)/abs((g)^(2/m) *
                                                              gss(x, clall[, 2], d)$wgss - (g + 1)^(2/m) *
                                                              gss(x, clall[, 3], d)$wgss)
      return(KL)
    }
    indice.ch <- function(x, cl, d = NULL, centrotypes = "centroids") {
      if (nb1.cl1 > 0) {
        CH <- NA
      }
      if (sum(c("centroids", "medoids") == centrotypes) ==
          0)
        stop("Wrong centrotypes argument")
      if ("medoids" == centrotypes && is.null(d))
        stop("For argument centrotypes = 'medoids' d cannot be null")
      if (!is.null(d)) {
        if (!is.matrix(d)) {
          d <- as.matrix(d)
        }
        row.names(d) <- row.names(x)
      }
      n <- length(cl)
      k <- max(cl)
      CH <- (gss(x, cl, d)$bgss/(k - 1))/(gss(x, cl, d)$wgss/(n -
                                                                k))
      return(CH)
    }
    indice.hart <- function(x, clall, d = NULL, centrotypes = "centroids") {
      if (sum(c("centroids", "medoids") == centrotypes) ==
          0)
        stop("Wrong centrotypes argument")
      if ("medoids" == centrotypes && is.null(d))
        stop("For argument centrotypes = 'medoids' d cannot be null")
      if (!is.null(d)) {
        if (!is.matrix(d)) {
          d <- as.matrix(d)
        }
        row.names(d) <- row.names(x)
      }
      n <- nrow(x)
      g <- max(clall[, 2])
      HART <- (gss(x, clall[, 2], d)$wgss/gss(x, clall[, 3], d)$wgss - 1) * (n - g - 1)
      return(HART)
    }
    indice.ball <- function(x, cl, d = NULL, centrotypes = "centroids") {
      wgssB <- gss(x, cl, d)$wgss
      qq <- max(cl)
      ball <- wgssB/qq
      return(ball)
    }
    indice.ratkowsky <- function(x, cl, d, centrotypes = "centroids") {
      qq <- max(cl)
      clsize <- table(cl)
      centers <- gss(x, cl, d)$centers
      varwithins <- varwithinss(x, centers, cl)
      zvargss <- vargss(x, clsize, varwithins)
      ratio <- mean(sqrt(zvargss$varbgss/zvargss$vartss))
      ratkowsky <- ratio/sqrt(qq)
      return(ratkowsky)
    }
    indice <- pmatch(index, c("kl", "ch", "hart", "ratkowsky",
                              "ball", "all"))
    if (is.na(indice))
      stop("invalid clustering index")
    if (indice == -1)
      stop("ambiguous index")
    vecallindex <- numeric(5)
    if (any(indice == 1) || (indice == 6))
      vecallindex[1] <- indice.kl(x, clall, d)
    if (any(indice == 2) || (indice == 6))
      vecallindex[2] <- indice.ch(x, cl = clall[, 2],
                                  d)
    if (any(indice == 3) || (indice == 6))
      vecallindex[3] <- indice.hart(x, clall, d)
    if (any(indice == 4) || (indice == 6))
      vecallindex[4] <- indice.ratkowsky(x, cl = cl1,
                                         d)
    if (any(indice == 5) || (indice == 6))
      vecallindex[5] <- indice.ball(x, cl = cl1, d)
    names(vecallindex) <- c("kl", "ch", "hart", "ratkowsky",
                            "ball")
    if (indice < 6)
      vecallindex <- vecallindex[indice]
    return(vecallindex)
  }
  Indice.cindex <- function(d, cl) {
    d <- data.matrix(d)
    DU <- 0
    r <- 0
    v_max <- array(1, max(cl))
    v_min <- array(1, max(cl))
    for (i in 1:max(cl)) {
      n <- sum(cl == i)
      if (n > 1) {
        t <- d[cl == i, cl == i]
        DU = DU + sum(t)/2
        v_max[i] = max(t)
        if (sum(t == 0) == n)
          v_min[i] <- min(t[t != 0])
        else v_min[i] <- 0
        r <- r + n * (n - 1)/2
      }
    }
    Dmin = min(v_min)
    Dmax = max(v_max)
    if (Dmin == Dmax)
      result <- NA
    else result <- (DU - r * Dmin)/(Dmax * r - Dmin * r)
    result
  }
  Indice.DB <- function(x, cl, d = NULL, centrotypes = "centroids",
                        p = 2, q = 2) {
    if (sum(c("centroids") == centrotypes) == 0)
      stop("Wrong centrotypes argument")
    if (!is.null(d)) {
      if (!is.matrix(d)) {
        d <- as.matrix(d)
      }
      row.names(d) <- row.names(x)
    }
    if (is.null(dim(x))) {
      dim(x) <- c(length(x), 1)
    }
    x <- as.matrix(x)
    n <- length(cl)
    k <- max(cl)
    dAm <- d
    centers <- matrix(nrow = k, ncol = ncol(x))
    if (centrotypes == "centroids") {
      for (i in 1:k) {
        for (j in 1:ncol(x)) {
          centers[i, j] <- mean(x[cl == i, j])
        }
      }
    }
    else {
      stop("wrong centrotypes argument")
    }
    S <- rep(0, k)
    for (i in 1:k) {
      ind <- (cl == i)
      if (sum(ind) > 1) {
        centerI <- centers[i, ]
        centerI <- rep(centerI, sum(ind))
        centerI <- matrix(centerI, nrow = sum(ind),
                          ncol = ncol(x), byrow = TRUE)
        S[i] <- mean(sqrt(apply((x[ind, ] - centerI)^2,
                                1, sum))^q)^(1/q)
      }
      else S[i] <- 0
    }
    M <- as.matrix(dist(centers, p = p))
    R <- array(Inf, c(k, k))
    r = rep(0, k)
    for (i in 1:k) {
      for (j in 1:k) {
        R[i, j] = (S[i] + S[j])/M[i, j]
      }
      r[i] = max(R[i, ][is.finite(R[i, ])])
    }
    DB = mean(r[is.finite(r)])
    resul <- list(DB = DB, r = r, R = R, d = M, S = S, centers = centers)
    resul
  }
  Indice.S <- function(d, cl) {
    d <- as.matrix(d)
    Si <- 0
    for (k in 1:max(cl)) {
      if ((sum(cl == k)) <= 1)
        Sil <- 1
      else {
        Sil <- 0
        for (i in 1:length(cl)) {
          if (cl[i] == k) {
            ai <- sum(d[i, cl == k])/(sum(cl == k) -
                                        1)
            dips <- NULL
            for (j in 1:max(cl)) if (cl[i] != j)
              if (sum(cl == j) != 1)
                dips <- cbind(dips, c((sum(d[i, cl ==
                                               j]))/(sum(cl == j))))
            else dips <- cbind(dips, c((sum(d[i, cl ==
                                                j]))))
            bi <- min(dips)
            Sil <- Sil + (bi - ai)/max(c(ai, bi))
          }
        }
      }
      Si <- Si + Sil
    }
    Si/length(cl)
  }
  Indice.Gap <- function(x, clall, reference.distribution = "unif",
                         B = 10, method = "ward.D2", d = NULL, centrotypes = "centroids") {
    GAP <- function(X, cl, referenceDistribution, B, method,
                    d, centrotypes) {
      set.seed(1)
      simgap <- function(Xvec) {
        ma <- max(Xvec)
        mi <- min(Xvec)
        set.seed(1)
        Xout <- runif(length(Xvec), min = mi, max = ma)
        return(Xout)
      }
      pcsim <- function(X, d, centrotypes) {
        if (centrotypes == "centroids") {
          Xmm <- apply(X, 2, mean)
        }
        for (k in (1:dim(X)[2])) {
          X[, k] <- X[, k] - Xmm[k]
        }
        ss <- svd(X)
        Xs <- X %*% ss$v
        Xnew <- apply(Xs, 2, simgap)
        Xt <- Xnew %*% t(ss$v)
        for (k in (1:dim(X)[2])) {
          Xt[, k] <- Xt[, k] + Xmm[k]
        }
        return(Xt)
      }
      if (is.null(dim(x))) {
        dim(x) <- c(length(x), 1)
      }
      ClassNr <- max(cl)
      Wk0 <- 0
      WkB <- matrix(0, 1, B)
      for (bb in (1:B)) {
        if (reference.distribution == "unif")
          Xnew <- apply(X, 2, simgap)
        else if (reference.distribution == "pc")
          Xnew <- pcsim(X, d, centrotypes)
        else stop("Wrong reference distribution type")
        if (bb == 1) {
          pp <- cl
          if (ClassNr == length(cl))
            pp2 <- 1:ClassNr
          else if (method == "k-means") {
            set.seed(1)
            pp2 <- kmeans(Xnew, ClassNr, 100)$cluster
          }
          else if (method == "single" || method == "complete" ||
                   method == "average" || method == "ward.D2" ||
                   method == "mcquitty" || method == "median" ||
                   method == "centroid" || method == "ward.D")
            pp2 <- cutree(hclust(dist(Xnew), method = method),
                          ClassNr)
          else stop("Wrong clustering method")
          if (ClassNr > 1) {
            for (zz in (1:ClassNr)) {
              Xuse <- X[pp == zz, ]
              Wk0 <- Wk0 + sum(diag(var(Xuse))) * (length(pp[pp ==
                                                               zz]) - 1)/(dim(X)[1] - ClassNr)
              Xuse2 <- Xnew[pp2 == zz, ]
              WkB[1, bb] <- WkB[1, bb] + sum(diag(var(Xuse2))) *
                (length(pp2[pp2 == zz]) - 1)/(dim(X)[1] -
                                                ClassNr)
            }
          }
          if (ClassNr == 1) {
            Wk0 <- sum(diag(var(X)))
            WkB[1, bb] <- sum(diag(var(Xnew)))
          }
        }
        if (bb > 1) {
          if (ClassNr == length(cl))
            pp2 <- 1:ClassNr
          else if (method == "k-means") {
            set.seed(1)
            pp2 <- kmeans(Xnew, ClassNr, 100)$cluster
          }
          else if (method == "single" || method == "complete" ||
                   method == "average" || method == "ward.D2" ||
                   method == "mcquitty" || method == "median" ||
                   method == "centroid" || method == "ward.D")
            pp2 <- cutree(hclust(dist(Xnew), method = method),
                          ClassNr)
          else stop("Wrong clustering method")
          if (ClassNr > 1) {
            for (zz in (1:ClassNr)) {
              Xuse2 <- Xnew[pp2 == zz, ]
              WkB[1, bb] <- WkB[1, bb] + sum(diag(var(Xuse2))) *
                length(pp2[pp2 == zz])/(dim(X)[1] -
                                          ClassNr)
            }
          }
          if (ClassNr == 1) {
            WkB[1, bb] <- sum(diag(var(Xnew)))
          }
        }
      }
      Sgap <- mean(log(WkB[1, ])) - log(Wk0)
      Sdgap <- sqrt(1 + 1/B) * sqrt(var(log(WkB[1, ]))) *
        sqrt((B - 1)/B)
      resul <- list(Sgap = Sgap, Sdgap = Sdgap)
      resul
    }
    if (sum(c("centroids", "medoids") == centrotypes) ==
        0)
      stop("Wrong centrotypes argument")
    if ("medoids" == centrotypes && is.null(d))
      stop("For argument centrotypes = 'medoids' d can not be null")
    if (!is.null(d)) {
      if (!is.matrix(d)) {
        d <- as.matrix(d)
      }
      row.names(d) <- row.names(x)
    }
    X <- as.matrix(x)
    gap1 <- GAP(X, clall[, 1], reference.distribution, B,
                method, d, centrotypes)
    gap <- gap1$Sgap
    gap2 <- GAP(X, clall[, 2], reference.distribution, B,
                method, d, centrotypes)
    diffu <- gap - (gap2$Sgap - gap2$Sdgap)
    resul <- list(gap = gap, diffu = diffu)
    resul
  }
  Index.sdindex <- function(x, clmax, cl) {
    x <- as.matrix(x)
    Alpha <- Dis(clmax, x)
    Scatt <- Average.scattering(cl, x)$scatt
    Dis0 <- Dis(cl, x)
    SD.indice <- Alpha * Scatt + Dis0
    return(SD.indice)
  }
  Index.SDbw <- function(x, cl) {
    x <- as.matrix(x)
    Scatt <- Average.scattering(cl, x)$scatt
    Dens.bw <- density.bw(cl, x)
    SDbw <- Scatt + Dens.bw
    return(SDbw)
  }
  Index.Dindex <- function(cl, x) {
    x <- as.matrix(x)
    distance <- density.clusters(cl, x)$distance
    n <- length(distance)
    S <- 0
    for (i in 1:n) S <- S + distance[i]
    inertieIntra <- S/n
    return(inertieIntra)
  }
  Index.dunn <- function(md, clusters, Data = NULL, method = "euclidean") {
    distance <- as.matrix(md)
    nc <- max(clusters)
    interClust <- matrix(NA, nc, nc)
    intraClust <- rep(NA, nc)
    for (i in 1:nc) {
      c1 <- which(clusters == i)
      for (j in i:nc) {
        if (j == i)
          intraClust[i] <- max(distance[c1, c1])
        if (j > i) {
          c2 <- which(clusters == j)
          interClust[i, j] <- min(distance[c1, c2])
        }
      }
    }
    dunn <- min(interClust, na.rm = TRUE)/max(intraClust)
    return(dunn)
  }
  for (nc in min_nc:max_nc) {
    if (any(method == 1) || (method == 2) || (method ==
                                              3) || (method == 4) || (method == 5) || (method ==
                                                                                       6) || (method == 7) || (method == 9)) {
      cl1 <- cutree(hc, k = nc)
      cl2 <- cutree(hc, k = nc + 1)
      clall <- cbind(cl1, cl2)
      clmax <- cutree(hc, k = max_nc)
      if (nc >= 2) {
        cl0 <- cutree(hc, k = nc - 1)
        clall1 <- cbind(cl0, cl1, cl2)
      }
      if (nc == 1) {
        cl0 <- rep(NA, nn)
        clall1 <- cbind(cl0, cl1, cl2)
      }
    }
    if (method == 8) {
      set.seed(1)
      cl2 <- kmeans(jeu, nc + 1, nstart=20)$cluster
      set.seed(1)
      clmax <- kmeans(jeu, max_nc, nstart=20)$cluster
      if (nc > 2) {
        set.seed(1)
        cl1 <- kmeans(jeu, nc, nstart=20)$cluster
        clall <- cbind(cl1, cl2)
        set.seed(1)
        cl0 <- kmeans(jeu, nc - 1, nstart=20)$cluster
        clall1 <- cbind(cl0, cl1, cl2)
      }
      if (nc == 2) {
        set.seed(1)
        cl1 <- kmeans(jeu, nc, nstart=20)$cluster
        clall <- cbind(cl1, cl2)
        cl0 <- rep(1, nn)
        clall1 <- cbind(cl0, cl1, cl2)
      }
      if (nc == 1) {
        stop("Number of clusters must be higher than 2")
      }
    }
    j <- table(cl1)
    s <- sum(j == 1)
    j2 <- table(cl2)
    s2 <- sum(j2 == 1)
    if (any(indice == 3) || (indice == 31) || (indice ==
                                               32)) {
      res[nc - min_nc + 1, 3] <- Indices.Traces(jeu, md,
                                                clall1, index = "hart")
    }
    if (any(indice == 4) || (indice == 31) || (indice ==
                                               32)) {
      res[nc - min_nc + 1, 4] <- Indices.WBT(x = jeu,
                                             cl = cl1, P = TT, s = ss, vv = vv)$ccc
    }
    if (any(indice == 5) || (indice == 31) || (indice ==
                                               32)) {
      res[nc - min_nc + 1, 5] <- Indices.WBT(x = jeu,
                                             cl = cl1, P = TT, s = ss, vv = vv)$scott
    }
    if (any(indice == 6) || (indice == 31) || (indice ==
                                               32)) {
      res[nc - min_nc + 1, 6] <- Indices.WBT(x = jeu,
                                             cl = cl1, P = TT, s = ss, vv = vv)$marriot
    }
    if (any(indice == 7) || (indice == 31) || (indice ==
                                               32)) {
      res[nc - min_nc + 1, 7] <- Indices.WBT(x = jeu,
                                             cl = cl1, P = TT, s = ss, vv = vv)$trcovw
    }
    if (any(indice == 8) || (indice == 31) || (indice ==
                                               32)) {
      res[nc - min_nc + 1, 8] <- Indices.WBT(x = jeu,
                                             cl = cl1, P = TT, s = ss, vv = vv)$tracew
    }
    if (any(indice == 9) || (indice == 31) || (indice ==
                                               32)) {
      res[nc - min_nc + 1, 9] <- Indices.WBT(x = jeu,
                                             cl = cl1, P = TT, s = ss, vv = vv)$friedman
    }
    if (any(indice == 10) || (indice == 31) || (indice ==
                                                32)) {
      res[nc - min_nc + 1, 10] <- Indices.WBT(x = jeu,
                                              cl = cl1, P = TT, s = ss, vv = vv)$rubin
    }
    if (any(indice == 14) || (indice == 31) || (indice ==
                                                32)) {
      res[nc - min_nc + 1, 14] <- Indices.WKWL(x = jeu,
                                               cl1 = cl1, cl2 = cl2)$duda
    }
    if (any(indice == 15) || (indice == 31) || (indice ==
                                                32)) {
      res[nc - min_nc + 1, 15] <- Indices.WKWL(x = jeu,
                                               cl1 = cl1, cl2 = cl2)$pseudot2
    }
    if (any(indice == 16) || (indice == 31) || (indice ==
                                                32)) {
      res[nc - min_nc + 1, 16] <- beale <- Indices.WKWL(x = jeu,
                                                        cl1 = cl1, cl2 = cl2)$beale
    }
    if (any(indice == 14) || (indice == 15) || (indice ==
                                                16) || (indice == 31) || (indice == 32)) {
      NM <- Indices.WKWL(x = jeu, cl1 = cl1, cl2 = cl2)$NM
      NK <- Indices.WKWL(x = jeu, cl1 = cl1, cl2 = cl2)$NK
      NL <- Indices.WKWL(x = jeu, cl1 = cl1, cl2 = cl2)$NL
      zz <- 3.2
      zzz <- zz * sqrt(2 * (1 - 8/((pi^2) * pp))/(NM *
                                                    pp))
      if (any(indice == 14) || (indice == 31) || (indice ==
                                                  32)) {
        resCritical[nc - min_nc + 1, 1] <- critValue <- 1 -
          (2/(pi * pp)) - zzz
      }
      if ((indice == 15) || (indice == 31) || (indice ==
                                               32)) {
        critValue <- 1 - (2/(pi * pp)) - zzz
        resCritical[nc - min_nc + 1, 2] <- ((1 - critValue)/critValue) *
          (NK + NL - 2)
      }
      if (any(indice == 16) || (indice == 31) || (indice ==
                                                  32)) {
        df2 <- (NM - 2) * pp
        resCritical[nc - min_nc + 1, 3] <- 1 - pf(beale,
                                                  pp, df2)
      }
    }
    if (any(indice == 18) || (indice == 31) || (indice ==
                                                32)) {
      res[nc - min_nc + 1, 18] <- Indices.Traces(jeu,
                                                 md, clall1, index = "ball")
    }
    if (any(indice == 19) || (indice == 31) || (indice ==
                                                32)) {
      res[nc - min_nc + 1, 19] <- Indice.ptbiserial(x = jeu,
                                                    md = md, cl1 = cl1)
    }
    if (any(indice == 20) || (indice == 32)) {
      if (method == 1) {
        resultSGAP <- Indice.Gap(x = jeu, clall = clall,
                                 reference.distribution = "unif", B = 10, method = "ward.D2",
                                 d = NULL, centrotypes = "centroids")
      }
      if (method == 2) {
        resultSGAP <- Indice.Gap(x = jeu, clall = clall,
                                 reference.distribution = "unif", B = 10, method = "single",
                                 d = NULL, centrotypes = "centroids")
      }
      if (method == 3) {
        resultSGAP <- Indice.Gap(x = jeu, clall = clall,
                                 reference.distribution = "unif", B = 10, method = "complete",
                                 d = NULL, centrotypes = "centroids")
      }
      if (method == 4) {
        resultSGAP <- Indice.Gap(x = jeu, clall = clall,
                                 reference.distribution = "unif", B = 10, method = "average",
                                 d = NULL, centrotypes = "centroids")
      }
      if (method == 5) {
        resultSGAP <- Indice.Gap(x = jeu, clall = clall,
                                 reference.distribution = "unif", B = 10, method = "mcquitty",
                                 d = NULL, centrotypes = "centroids")
      }
      if (method == 6) {
        resultSGAP <- Indice.Gap(x = jeu, clall = clall,
                                 reference.distribution = "unif", B = 10, method = "median",
                                 d = NULL, centrotypes = "centroids")
      }
      if (method == 7) {
        resultSGAP <- Indice.Gap(x = jeu, clall = clall,
                                 reference.distribution = "unif", B = 10, method = "centroid",
                                 d = NULL, centrotypes = "centroids")
      }
      if (method == 9) {
        resultSGAP <- Indice.Gap(x = jeu, clall = clall,
                                 reference.distribution = "unif", B = 10, method = "ward.D",
                                 d = NULL, centrotypes = "centroids")
      }
      if (method == 8) {
        resultSGAP <- Indice.Gap(x = jeu, clall = clall,
                                 reference.distribution = "unif", B = 10, method = "k-means",
                                 d = NULL, centrotypes = "centroids")
      }
      res[nc - min_nc + 1, 20] <- resultSGAP$gap
      resCritical[nc - min_nc + 1, 4] <- resultSGAP$diffu
    }
    if (nc >= 2) {
      if (any(indice == 1) || (indice == 31) || (indice ==
                                                 32)) {
        res[nc - min_nc + 1, 1] <- Indices.Traces(jeu,
                                                  md, clall1, index = "kl")
      }
      if (any(indice == 2) || (indice == 31) || (indice ==
                                                 32)) {
        res[nc - min_nc + 1, 2] <- Indices.Traces(jeu,
                                                  md, clall1, index = "ch")
      }
      if (any(indice == 11) || (indice == 31) || (indice ==
                                                  32)) {
        res[nc - min_nc + 1, 11] <- Indice.cindex(d = md,
                                                  cl = cl1)
      }
      if (any(indice == 12) || (indice == 31) || (indice ==
                                                  32)) {
        res[nc - min_nc + 1, 12] <- Indice.DB(x = jeu,
                                              cl = cl1, d = NULL, centrotypes = "centroids",
                                              p = 2, q = 2)$DB
      }
      if (any(indice == 13) || (indice == 31) || (indice ==
                                                  32)) {
        res[nc - min_nc + 1, 13] <- Indice.S(d = md,
                                             cl = cl1)
      }
      if (any(indice == 17) || (indice == 31) || (indice ==
                                                  32)) {
        res[nc - min_nc + 1, 17] <- Indices.Traces(jeu,
                                                   md, clall1, index = "ratkowsky")
      }
      if (any(indice == 21) || (indice == 31) || (indice ==
                                                  32)) {
        res[nc - min_nc + 1, 21] <- Index.15and28(cl1 = cl1,
                                                  cl2 = cl2, md = md)$frey
      }
      if (any(indice == 22) || (indice == 31) || (indice ==
                                                  32)) {
        res[nc - min_nc + 1, 22] <- Index.15and28(cl1 = cl1,
                                                  cl2 = cl2, md = md)$mcclain
      }
      if (any(indice == 23) || (indice == 32)) {
        res[nc - min_nc + 1, 23] <- Index.sPlussMoins(cl1 = cl1,
                                                      md = md)$gamma
      }
      if (any(indice == 24) || (indice == 32)) {
        res[nc - min_nc + 1, 24] <- Index.sPlussMoins(cl1 = cl1,
                                                      md = md)$gplus
      }
      if (any(indice == 25) || (indice == 32)) {
        res[nc - min_nc + 1, 25] <- Index.sPlussMoins(cl1 = cl1,
                                                      md = md)$tau
      }
      if (any(indice == 26) || (indice == 31) || (indice ==
                                                  32)) {
        res[nc - min_nc + 1, 26] <- Index.dunn(md, cl1,
                                               Data = jeu, method = NULL)
      }
      if (any(indice == 27) || (indice == 31) || (indice ==
                                                  32)) {
        res[nc - min_nc + 1, 27] <- Index.Hubert(jeu,
                                                 cl1)
      }
      if (any(indice == 28) || (indice == 31) || (indice ==
                                                  32)) {
        res[nc - min_nc + 1, 28] <- Index.sdindex(jeu,
                                                  clmax, cl1)
      }
      if (any(indice == 29) || (indice == 31) || (indice ==
                                                  32)) {
        res[nc - min_nc + 1, 29] <- Index.Dindex(cl1,
                                                 jeu)
      }
      if (any(indice == 30) || (indice == 31) || (indice ==
                                                  32)) {
        res[nc - min_nc + 1, 30] <- Index.SDbw(jeu,
                                               cl1)
      }
    }
    else {
      res[nc - min_nc + 1, 1] <- NA
      res[nc - min_nc + 1, 2] <- NA
      res[nc - min_nc + 1, 11] <- NA
      res[nc - min_nc + 1, 12] <- NA
      res[nc - min_nc + 1, 13] <- NA
      res[nc - min_nc + 1, 17] <- NA
      res[nc - min_nc + 1, 21] <- NA
      res[nc - min_nc + 1, 22] <- NA
      res[nc - min_nc + 1, 23] <- NA
      res[nc - min_nc + 1, 24] <- NA
      res[nc - min_nc + 1, 25] <- NA
      res[nc - min_nc + 1, 26] <- NA
      res[nc - min_nc + 1, 27] <- NA
      res[nc - min_nc + 1, 28] <- NA
      res[nc - min_nc + 1, 29] <- NA
      res[nc - min_nc + 1, 30] <- NA
    }
  }
  nc.KL <- indice.KL <- 0
  if (any(indice == 1) || (indice == 31) || (indice == 32)) {
    nc.KL <- (min_nc:max_nc)[which.max(res[, 1])]
    indice.KL <- max(res[, 1], na.rm = TRUE)
    best.nc <- nc.KL
  }
  nc.CH <- indice.CH <- 0
  if (any(indice == 2) || (indice == 31) || (indice == 32)) {
    nc.CH <- (min_nc:max_nc)[which.max(res[, 2])]
    indice.CH <- max(res[, 2], na.rm = TRUE)
    best.nc <- nc.CH
  }
  nc.CCC <- indice.CCC <- 0
  if (any(indice == 4) || (indice == 31) || (indice == 32)) {
    nc.CCC <- (min_nc:max_nc)[which.max(res[, 4])]
    indice.CCC <- max(res[, 4], na.rm = TRUE)
    best.nc <- nc.CCC
  }
  nc.DB <- indice.DB <- 0
  if (any(indice == 12) || (indice == 31) || (indice == 32)) {
    nc.DB <- (min_nc:max_nc)[which.min(res[, 12])]
    indice.DB <- min(res[, 12], na.rm = TRUE)
    best.nc <- nc.DB
  }
  nc.Silhouette <- indice.Silhouette <- 0
  if (any(indice == 13) || (indice == 31) || (indice == 32)) {
    nc.Silhouette <- (min_nc:max_nc)[which.max(res[, 13])]
    indice.Silhouette <- max(res[, 13], na.rm = TRUE)
    best.nc <- nc.Silhouette
  }
  nc.Gap <- indice.Gap <- 0
  if (any(indice == 20) || (indice == 32)) {
    found <- FALSE
    for (ncG in min_nc:max_nc) {
      if ((resCritical[ncG - min_nc + 1, 4] >= 0) && (!found)) {
        ncGap <- ncG
        indiceGap <- res[ncG - min_nc + 1, 20]
        found <- TRUE
      }
    }
    if (found) {
      nc.Gap <- ncGap
      indice.Gap <- indiceGap
      best.nc <- nc.Gap
    }
    else {
      nc.Gap <- NA
      indice.Gap <- NA
    }
  }
  nc.Duda <- indice.Duda <- 0
  if (any(indice == 14) || (indice == 31) || (indice == 32)) {
    foundDuda <- FALSE
    for (ncD in min_nc:max_nc) {
      if ((res[ncD - min_nc + 1, 14] >= resCritical[ncD -
                                                    min_nc + 1, 1]) && (!foundDuda)) {
        ncDuda <- ncD
        indiceDuda <- res[ncD - min_nc + 1, 14]
        foundDuda <- TRUE
      }
    }
    if (foundDuda) {
      nc.Duda <- ncDuda
      indice.Duda <- indiceDuda
      best.nc <- nc.Duda
    }
    else {
      nc.Duda <- NA
      indice.Duda <- NA
    }
  }
  nc.Pseudo <- indice.Pseudo <- 0
  if (any(indice == 15) || (indice == 31) || (indice == 32)) {
    foundPseudo <- FALSE
    for (ncP in min_nc:max_nc) {
      if ((res[ncP - min_nc + 1, 15] <= resCritical[ncP -
                                                    min_nc + 1, 2]) && (!foundPseudo)) {
        ncPseudo <- ncP
        indicePseudo <- res[ncP - min_nc + 1, 15]
        foundPseudo <- TRUE
      }
    }
    if (foundPseudo) {
      nc.Pseudo <- ncPseudo
      indice.Pseudo <- indicePseudo
      best.nc <- nc.Pseudo
    }
    else {
      nc.Pseudo <- NA
      indice.Pseudo <- NA
    }
  }
  nc.Beale <- indice.Beale <- 0
  if (any(indice == 16) || (indice == 31) || (indice == 32)) {
    foundBeale <- FALSE
    for (ncB in min_nc:max_nc) {
      if ((resCritical[ncB - min_nc + 1, 3] >= alphaBeale) &&
          (!foundBeale)) {
        ncBeale <- ncB
        indiceBeale <- res[ncB - min_nc + 1, 16]
        foundBeale <- TRUE
      }
    }
    if (foundBeale) {
      nc.Beale <- ncBeale
      indice.Beale <- indiceBeale
      best.nc <- nc.Beale
    }
    else {
      nc.Beale <- NA
      indice.Beale <- NA
    }
  }
  nc.ptbiserial <- indice.ptbiserial <- 0
  if (any(indice == 19) || (indice == 31) || (indice == 32)) {
    nc.ptbiserial <- (min_nc:max_nc)[which.max(res[, 19])]
    indice.ptbiserial <- max(res[, 19], na.rm = TRUE)
    best.nc <- nc.ptbiserial
  }
  foundNC <- foundIndice <- numeric(0)
  nc.Frey <- indice.Frey <- 0
  if (any(indice == 21) || (indice == 31) || (indice == 32)) {
    foundFrey <- FALSE
    i <- 1
    for (ncF in min_nc:max_nc) {
      if (res[ncF - min_nc + 1, 21] < 1) {
        ncFrey <- ncF - 1
        indiceFrey <- res[ncF - 1 - min_nc + 1, 21]
        foundFrey <- TRUE
        foundNC[i] <- ncFrey
        foundIndice[i] <- indiceFrey
        i <- i + 1
      }
    }
    if (foundFrey) {
      nc.Frey <- foundNC[1]
      indice.Frey <- foundIndice[1]
      best.nc <- nc.Frey
    }
    else {
      nc.Frey <- NA
      indice.Frey <- NA
      print(paste("Frey index : No clustering structure in this data set"))
    }
  }
  nc.McClain <- indice.McClain <- 0
  if (any(indice == 22) || (indice == 31) || (indice == 32)) {
    nc.McClain <- (min_nc:max_nc)[which.min(res[, 22])]
    indice.McClain <- min(res[, 22], na.rm = TRUE)
    best.nc <- nc.McClain
  }
  nc.Gamma <- indice.Gamma <- 0
  if (any(indice == 23) || (indice == 31) || (indice == 32)) {
    nc.Gamma <- (min_nc:max_nc)[which.max(res[, 23])]
    indice.Gamma <- max(res[, 23], na.rm = TRUE)
    best.nc <- nc.Gamma
  }
  nc.Gplus <- indice.Gplus <- 0
  if (any(indice == 24) || (indice == 31) || (indice == 32)) {
    nc.Gplus <- (min_nc:max_nc)[which.min(res[, 24])]
    indice.Gplus <- min(res[, 24], na.rm = TRUE)
    best.nc <- nc.Gplus
  }
  nc.Tau <- indice.Tau <- 0
  if (any(indice == 25) || (indice == 31) || (indice == 32)) {
    nc.Tau <- (min_nc:max_nc)[which.max(res[, 25])]
    indice.Tau <- max(res[, 25], na.rm = TRUE)
    best.nc <- nc.Tau
  }
  if ((indice == 3) || (indice == 5) || (indice == 6) || (indice ==
                                                          7) || (indice == 8) || (indice == 9) || (indice == 10) ||
      (indice == 18) || (indice == 27) || (indice == 11) ||
      (indice == 29) || (indice == 31) || (indice == 32)) {
    DiffLev <- array(0, c(max_nc - min_nc + 1, 12))
    DiffLev[, 1] <- min_nc:max_nc
    for (nc3 in min_nc:max_nc) {
      if (nc3 == min_nc) {
        DiffLev[nc3 - min_nc + 1, 2] <- abs(res[nc3 -
                                                  min_nc + 1, 3] - NA)
        DiffLev[nc3 - min_nc + 1, 3] <- abs(res[nc3 -
                                                  min_nc + 1, 5] - NA)
        DiffLev[nc3 - min_nc + 1, 4] <- abs(res[nc3 -
                                                  min_nc + 1, 6] - NA)
        DiffLev[nc3 - min_nc + 1, 5] <- abs(res[nc3 -
                                                  min_nc + 1, 7] - NA)
        DiffLev[nc3 - min_nc + 1, 6] <- abs(res[nc3 -
                                                  min_nc + 1, 8] - NA)
        DiffLev[nc3 - min_nc + 1, 7] <- abs(res[nc3 -
                                                  min_nc + 1, 9] - NA)
        DiffLev[nc3 - min_nc + 1, 8] <- abs(res[nc3 -
                                                  min_nc + 1, 10] - NA)
        DiffLev[nc3 - min_nc + 1, 9] <- abs(res[nc3 -
                                                  min_nc + 1, 18] - NA)
        DiffLev[nc3 - min_nc + 1, 10] <- abs(res[nc3 -
                                                   min_nc + 1, 27] - NA)
        DiffLev[nc3 - min_nc + 1, 12] <- abs(res[nc3 -
                                                   min_nc + 1, 29] - NA)
      }
      else {
        if (nc3 == max_nc) {
          DiffLev[nc3 - min_nc + 1, 2] <- abs(res[nc3 -
                                                    min_nc + 1, 3] - res[nc3 - min_nc, 3])
          DiffLev[nc3 - min_nc + 1, 3] <- abs(res[nc3 -
                                                    min_nc + 1, 5] - res[nc3 - min_nc, 5])
          DiffLev[nc3 - min_nc + 1, 4] <- abs(res[nc3 -
                                                    min_nc + 1, 6] - NA)
          DiffLev[nc3 - min_nc + 1, 5] <- abs(res[nc3 -
                                                    min_nc + 1, 7] - res[nc3 - min_nc, 7])
          DiffLev[nc3 - min_nc + 1, 6] <- abs(res[nc3 -
                                                    min_nc + 1, 8] - NA)
          DiffLev[nc3 - min_nc + 1, 7] <- abs(res[nc3 -
                                                    min_nc + 1, 9] - res[nc3 - min_nc, 9])
          DiffLev[nc3 - min_nc + 1, 8] <- abs(res[nc3 -
                                                    min_nc + 1, 10] - NA)
          DiffLev[nc3 - min_nc + 1, 9] <- abs(res[nc3 -
                                                    min_nc + 1, 18] - res[nc3 - min_nc, 18])
          DiffLev[nc3 - min_nc + 1, 10] <- abs(res[nc3 -
                                                     min_nc + 1, 27] - NA)
          DiffLev[nc3 - min_nc + 1, 12] <- abs(res[nc3 -
                                                     min_nc + 1, 29] - NA)
        }
        else {
          DiffLev[nc3 - min_nc + 1, 2] <- abs(res[nc3 -
                                                    min_nc + 1, 3] - res[nc3 - min_nc, 3])
          DiffLev[nc3 - min_nc + 1, 3] <- abs(res[nc3 -
                                                    min_nc + 1, 5] - res[nc3 - min_nc, 5])
          DiffLev[nc3 - min_nc + 1, 4] <- ((res[nc3 -
                                                  min_nc + 2, 6] - res[nc3 - min_nc + 1, 6]) -
                                             (res[nc3 - min_nc + 1, 6] - res[nc3 - min_nc,
                                                                             6]))
          DiffLev[nc3 - min_nc + 1, 5] <- abs(res[nc3 -
                                                    min_nc + 1, 7] - res[nc3 - min_nc, 7])
          DiffLev[nc3 - min_nc + 1, 6] <- ((res[nc3 -
                                                  min_nc + 2, 8] - res[nc3 - min_nc + 1, 8]) -
                                             (res[nc3 - min_nc + 1, 8] - res[nc3 - min_nc,
                                                                             8]))
          DiffLev[nc3 - min_nc + 1, 7] <- abs(res[nc3 -
                                                    min_nc + 1, 9] - res[nc3 - min_nc, 9])
          DiffLev[nc3 - min_nc + 1, 8] <- ((res[nc3 -
                                                  min_nc + 2, 10] - res[nc3 - min_nc + 1,
                                                                        10]) - (res[nc3 - min_nc + 1, 10] - res[nc3 -
                                                                                                                  min_nc, 10]))
          DiffLev[nc3 - min_nc + 1, 9] <- abs(res[nc3 -
                                                    min_nc + 1, 18] - res[nc3 - min_nc, 18])
          DiffLev[nc3 - min_nc + 1, 10] <- abs((res[nc3 -
                                                      min_nc + 1, 27] - res[nc3 - min_nc, 27]))
          DiffLev[nc3 - min_nc + 1, 12] <- ((res[nc3 -
                                                   min_nc + 2, 29] - res[nc3 - min_nc + 1,
                                                                         29]) - (res[nc3 - min_nc + 1, 29] - res[nc3 -
                                                                                                                   min_nc, 29]))
        }
      }
    }
  }
  nc.Hartigan <- indice.Hartigan <- 0
  if (any(indice == 3) || (indice == 31) || (indice == 32)) {
    nc.Hartigan <- DiffLev[, 1][which.max(DiffLev[, 2])]
    indice.Hartigan <- max(DiffLev[, 2], na.rm = TRUE)
    best.nc <- nc.Hartigan
  }
  nc.Ratkowsky <- indice.Ratkowsky <- 0
  if (any(indice == 17) || (indice == 31) || (indice == 32)) {
    nc.Ratkowsky <- (min_nc:max_nc)[which.max(res[, 17])]
    indice.Ratkowsky <- max(res[, 17], na.rm = TRUE)
    best.nc <- nc.Ratkowsky
  }
  nc.cindex <- indice.cindex <- 0
  if (any(indice == 11) || (indice == 31) || (indice == 32)) {
    nc.cindex <- (min_nc:max_nc)[which.min(res[, 11])]
    indice.cindex <- min(res[, 11], na.rm = TRUE)
    best.nc <- nc.cindex
  }
  nc.Scott <- indice.Scott <- 0
  if (any(indice == 5) || (indice == 31) || (indice == 32)) {
    nc.Scott <- DiffLev[, 1][which.max(DiffLev[, 3])]
    indice.Scott <- max(DiffLev[, 3], na.rm = TRUE)
    best.nc <- nc.Scott
  }
  nc.Marriot <- indice.Marriot <- 0
  if (any(indice == 6) || (indice == 31) || (indice == 32)) {
    nc.Marriot <- DiffLev[, 1][which.max(DiffLev[, 4])]
    round(nc.Marriot, digits = 1)
    indice.Marriot <- max(DiffLev[, 4], na.rm = TRUE)
    best.nc <- nc.Marriot
  }
  nc.TrCovW <- indice.TrCovW <- 0
  if (any(indice == 7) || (indice == 31) || (indice == 32)) {
    nc.TrCovW <- DiffLev[, 1][which.max(DiffLev[, 5])]
    indice.TrCovW <- max(DiffLev[, 5], na.rm = TRUE)
    best.nc <- nc.TrCovW
  }
  nc.TraceW <- indice.TraceW <- 0
  if (any(indice == 8) || (indice == 31) || (indice == 32)) {
    nc.TraceW <- DiffLev[, 1][which.max(DiffLev[, 6])]
    indice.TraceW <- max(DiffLev[, 6], na.rm = TRUE)
    best.nc <- nc.TraceW
  }
  nc.Friedman <- indice.Friedman <- 0
  if (any(indice == 9) || (indice == 31) || (indice == 32)) {
    nc.Friedman <- DiffLev[, 1][which.max(DiffLev[, 7])]
    indice.Friedman <- max(DiffLev[, 7], na.rm = TRUE)
    best.nc <- nc.Friedman
  }
  nc.Rubin <- indice.Rubin <- 0
  if (any(indice == 10) || (indice == 31) || (indice == 32)) {
    nc.Rubin <- DiffLev[, 1][which.min(DiffLev[, 8])]
    indice.Rubin <- min(DiffLev[, 8], na.rm = TRUE)
    best.nc <- nc.Rubin
  }
  nc.Ball <- indice.Ball <- 0
  if (any(indice == 18) || (indice == 31) || (indice == 32)) {
    nc.Ball <- DiffLev[, 1][which.max(DiffLev[, 9])]
    indice.Ball <- max(DiffLev[, 9], na.rm = TRUE)
    best.nc <- nc.Ball
  }
  nc.Dunn <- indice.Dunn <- 0
  if (any(indice == 26) || (indice == 31) || (indice == 32)) {
    nc.Dunn <- (min_nc:max_nc)[which.max(res[, 26])]
    indice.Dunn <- max(res[, 26], na.rm = TRUE)
    best.nc <- nc.Dunn
  }
  nc.Hubert <- indice.Hubert <- 0
  if (any(indice == 27) || (indice == 31) || (indice == 32)) {
    nc.Hubert <- 0
    indice.Hubert <- 0
    par(mfrow = c(1, 2))
    plot(x_axis, res[, 27], tck = 0, type = "b", col = "red",
         xlab = expression(paste("Number of clusters ")),
         ylab = expression(paste("Hubert Statistic values")))
    plot(DiffLev[, 1], DiffLev[, 10], tck = 0, type = "b",
         col = "blue", xlab = expression(paste("Number of clusters ")),
         ylab = expression(paste("Hubert statistic second differences")))
    cat(paste("*** : The Hubert index is a graphical method of determining the number of clusters.\n                In the plot of Hubert index, we seek a significant knee that corresponds to a \n                significant increase of the value of the measure i.e the significant peak in Hubert\n                index second differences plot.",
              "\n", "\n"))
  }
  nc.sdindex <- indice.sdindex <- 0
  if (any(indice == 28) || (indice == 31) || (indice == 32)) {
    nc.sdindex <- (min_nc:max_nc)[which.min(res[, 28])]
    indice.sdindex <- min(res[, 28], na.rm = TRUE)
    best.nc <- nc.sdindex
  }
  nc.Dindex <- indice.Dindex <- 0
  if (any(indice == 29) || (indice == 31) || (indice == 32)) {
    nc.Dindex <- 0
    indice.Dindex <- 0
    par(mfrow = c(1, 2))
    plot(x_axis, res[, 29], tck = 0, type = "b", col = "red",
         xlab = expression(paste("Number of clusters ")),
         ylab = expression(paste("Dindex Values")))
    plot(DiffLev[, 1], DiffLev[, 12], tck = 0, type = "b",
         col = "blue", xlab = expression(paste("Number of clusters ")),
         ylab = expression(paste("Second differences Dindex Values")))
    cat(paste("*** : The D index is a graphical method of determining the number of clusters. \n                In the plot of D index, we seek a significant knee (the significant peak in Dindex\n                second differences plot) that corresponds to a significant increase of the value of\n                the measure.",
              "\n", "\n"))
  }
  nc.SDbw <- indice.SDbw <- 0
  if (any(indice == 30) || (indice == 31) || (indice == 32)) {
    nc.SDbw <- (min_nc:max_nc)[which.min(res[, 30])]
    indice.SDbw <- min(res[, 30], na.rm = TRUE)
    best.nc <- nc.SDbw
  }
  if (indice < 31) {
    res <- res[, c(indice)]
    if (indice == 14) {
      resCritical <- resCritical[, 1]
    }
    if (indice == 15) {
      resCritical <- resCritical[, 2]
    }
    if (indice == 16) {
      resCritical <- resCritical[, 3]
    }
    if (indice == 20) {
      resCritical <- resCritical[, 4]
    }
  }
  if (indice == 31) {
    res <- res[, c(1:19, 21:22, 26:30)]
    resCritical <- resCritical[, c(1:3)]
  }
  if (any(indice == 20) || (indice == 23) || (indice == 24) ||
      (indice == 25) || (indice == 32)) {
    results <- c(nc.KL, indice.KL, nc.CH, indice.CH, nc.Hartigan,
                 indice.Hartigan, nc.CCC, indice.CCC, nc.Scott, indice.Scott,
                 nc.Marriot, indice.Marriot, nc.TrCovW, indice.TrCovW,
                 nc.TraceW, indice.TraceW, nc.Friedman, indice.Friedman,
                 nc.Rubin, indice.Rubin, nc.cindex, indice.cindex,
                 nc.DB, indice.DB, nc.Silhouette, indice.Silhouette,
                 nc.Duda, indice.Duda, nc.Pseudo, indice.Pseudo,
                 nc.Beale, indice.Beale, nc.Ratkowsky, indice.Ratkowsky,
                 nc.Ball, indice.Ball, nc.ptbiserial, indice.ptbiserial,
                 nc.Gap, indice.Gap, nc.Frey, indice.Frey, nc.McClain,
                 indice.McClain, nc.Gamma, indice.Gamma, nc.Gplus,
                 indice.Gplus, nc.Tau, indice.Tau, nc.Dunn, indice.Dunn,
                 nc.Hubert, indice.Hubert, nc.sdindex, indice.sdindex,
                 nc.Dindex, indice.Dindex, nc.SDbw, indice.SDbw)
    results1 <- matrix(c(results), nrow = 2, ncol = 30)
    resultats <- matrix(c(results), nrow = 2, ncol = 30,
                        dimnames = list(c("Number_clusters", "Value_Index"),
                                        c("KL", "CH", "Hartigan", "CCC", "Scott", "Marriot",
                                          "TrCovW", "TraceW", "Friedman", "Rubin", "Cindex",
                                          "DB", "Silhouette", "Duda", "PseudoT2", "Beale",
                                          "Ratkowsky", "Ball", "PtBiserial", "Gap",
                                          "Frey", "McClain", "Gamma", "Gplus", "Tau",
                                          "Dunn", "Hubert", "SDindex", "Dindex", "SDbw")))
  }
  else {
    results <- c(nc.KL, indice.KL, nc.CH, indice.CH, nc.Hartigan,
                 indice.Hartigan, nc.CCC, indice.CCC, nc.Scott, indice.Scott,
                 nc.Marriot, indice.Marriot, nc.TrCovW, indice.TrCovW,
                 nc.TraceW, indice.TraceW, nc.Friedman, indice.Friedman,
                 nc.Rubin, indice.Rubin, nc.cindex, indice.cindex,
                 nc.DB, indice.DB, nc.Silhouette, indice.Silhouette,
                 nc.Duda, indice.Duda, nc.Pseudo, indice.Pseudo,
                 nc.Beale, indice.Beale, nc.Ratkowsky, indice.Ratkowsky,
                 nc.Ball, indice.Ball, nc.ptbiserial, indice.ptbiserial,
                 nc.Frey, indice.Frey, nc.McClain, indice.McClain,
                 nc.Dunn, indice.Dunn, nc.Hubert, indice.Hubert,
                 nc.sdindex, indice.sdindex, nc.Dindex, indice.Dindex,
                 nc.SDbw, indice.SDbw)
    results1 <- matrix(c(results), nrow = 2, ncol = 26)
    resultats <- matrix(c(results), nrow = 2, ncol = 26,
                        dimnames = list(c("Number_clusters", "Value_Index"),
                                        c("KL", "CH", "Hartigan", "CCC", "Scott", "Marriot",
                                          "TrCovW", "TraceW", "Friedman", "Rubin", "Cindex",
                                          "DB", "Silhouette", "Duda", "PseudoT2", "Beale",
                                          "Ratkowsky", "Ball", "PtBiserial", "Frey",
                                          "McClain", "Dunn", "Hubert", "SDindex", "Dindex",
                                          "SDbw")))
  }
  if (any(indice <= 20) || (indice == 23) || (indice == 24) ||
      (indice == 25)) {
    resultats <- resultats[, c(indice)]
  }
  if (any(indice == 21) || (indice == 22)) {
    indice3 <- indice - 1
    resultats <- resultats[, c(indice3)]
  }
  if (any(indice == 26) || (indice == 27) || (indice == 28) ||
      (indice == 29) || (indice == 30)) {
    indice4 <- indice - 4
    resultats <- resultats[, c(indice4)]
  }
  resultats <- round(resultats, digits = 4)
  res <- round(res, digits = 4)
  resCritical <- round(resCritical, digits = 4)
  if (any(indice == 31) || (indice == 32)) {
    cat("*******************************************************************",
        "\n")
    cat("* Among all indices:                                               ",
        "\n")
    BestCluster <- results1[1, ]
    c = 0
    for (i in min.nc:max.nc) {
      vect <- which(BestCluster == i)
      if (length(vect) > 0)
        cat("*", length(vect), "proposed", i, "as the best number of clusters",
            "\n")
      if (c < length(vect)) {
        j = i
        c <- length(vect)
      }
    }
    cat("\n", "                  ***** Conclusion *****                           ",
        "\n", "\n")
    cat("* According to the majority rule, the best number of clusters is ",
        j, "\n", "\n", "\n")
    cat("*******************************************************************",
        "\n")
    if (any(method == 1) || (method == 2) || (method ==
                                              3) || (method == 4) || (method == 5) || (method ==
                                                                                       6) || (method == 7) || (method == 9))
      partition <- cutree(hc, k = j)
    else {
      set.seed(1)
      partition <- kmeans(jeu, j)$cluster
    }
  }
  if (any(indice == 1) || (indice == 2) || (indice == 3) ||
      (indice == 4) || (indice == 5) || (indice == 6) || (indice ==
                                                          7) || (indice == 8) || (indice == 9) || (indice == 10) ||
      (indice == 11) || (indice == 12) || (indice == 13) ||
      (indice == 14) || (indice == 15) || (indice == 16) ||
      (indice == 17) || (indice == 18) || (indice == 19) ||
      (indice == 20) || (indice == 21) || (indice == 22) ||
      (indice == 23) || (indice == 24) || (indice == 25) ||
      (indice == 26) || (indice == 28) || (indice == 30)) {
    if (any(method == 1) || (method == 2) || (method ==
                                              3) || (method == 4) || (method == 5) || (method ==
                                                                                       6) || (method == 7) || (method == 9))
      partition <- cutree(hc, k = best.nc)
    else {
      set.seed(1)
      partition <- kmeans(jeu, best.nc)$cluster
    }
  }
  if ((indice == 14) || (indice == 15) || (indice == 16) ||
      (indice == 20) || (indice == 31) || (indice == 32)) {
    results.final <- list(All.index = res, All.CriticalValues = resCritical,
                          Best.nc = resultats, Best.partition = partition)
  }
  if ((indice == 27) || (indice == 29))
    results.final <- list(All.index = res)
  if (any(indice == 1) || (indice == 2) || (indice == 3) ||
      (indice == 4) || (indice == 5) || (indice == 6) || (indice ==
                                                          7) || (indice == 8) || (indice == 9) || (indice == 10) ||
      (indice == 11) || (indice == 12) || (indice == 13) ||
      (indice == 17) || (indice == 18) || (indice == 19) ||
      (indice == 21) || (indice == 22) || (indice == 23) ||
      (indice == 24) || (indice == 25) || (indice == 26) ||
      (indice == 28) || (indice == 30))
    results.final <- list(All.index = res, Best.nc = resultats,
                          Best.partition = partition)
  return(results.final)
}


resample_sensitivity<-function(x, K.try = 2:10,n.resample,cut.off=0.8){
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
  if(max(result_data)<cut.off){
    maxk<-1
  }else{
    maxk<-K.try[which.max(result_data)]
  }

  return(maxk)
}
taylor.sparse_two_supplement = function(x, wbounds_list,trim=0.05,k_vector=c(2,3,4),n.resample=25){
  n = nrow(x)
  sub.n = round(n*0.7)
  #M.try = exp(seq(log(1.2), log(sqrt(ncol(x)) * 0.9),
  #               len = m.try))###???
  score_wbounds<-function(wbounds,x,k){
    print(wbounds)
    result.subsampling = lapply(1:n.resample, function(b){
      index.subsample = sample(1:nrow(x), sub.n, replace = F)
      xb = x[index.subsample,]
      #km.out <- kmeans(xb, centers = K, nstart = 100)
      #hc = hclust(dist(xb), method = "complete", members = NULL)
      #cl<-cmeans(xb,3,m=K)
      #set.seed(b)
      b = KMeansSparseCluster1(xb, K=k, wbounds = wbounds, nstart = 20, ###################
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
    o.result <- KMeansSparseCluster1(x, K=k, wbounds = wbounds, nstart = 100, ##############
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

    # stat.calc = function(i){
    #   if(all(c(1,0)%in%o.mtx[i,i:nr])){
    #     (mean(r.mtx[i,i:nr][which(o.mtx[i,i:nr] == 0)]) +
    #        mean(r.mtx[i,i:nr][which(o.mtx[i,i:nr] == 1)])) - 1
    #  }
    #   else if(1%in%o.mtx[i,i:nr]){
    #     mean(r.mtx[i,i:nr][which(o.mtx[i,i:nr] == 1)])
    #   }
    #   else {
    #     mean(r.mtx[i,i:nr][which(o.mtx[i,i:nr] == 0)])
    #   }
    # }
    # stat = rev(sapply(1:nr, stat.calc))
    # barplot(stat)

    ##so the original order is increasing
    stat.calc = function(i){
      if(all(c(1,0)%in%o.mtx[i,i:nr])){
        (mean(r.mtx[i,i:nr][which(o.mtx[i,i:nr] == 0)]) +
           mean(r.mtx[i,i:nr][which(o.mtx[i,i:nr] == 1)])) - 1
      }
      else if(1%in%o.mtx[i,i:nr]){
        mean(r.mtx[i,i:nr][which(o.mtx[i,i:nr] == 1)])
      }
      else {
        mean(r.mtx[i,i:nr][which(o.mtx[i,i:nr] == 0)])
      }
    }
    stat = rev(sapply(1:nr, stat.calc))
    barplot(stat)
    per<-1-trim
    clus_score = mean(stat[1:round(nr*per)])


    feature = sapply(result.subsampling, function(x) x$feature, simplify = "matrix")

    f.mtx  =  apply(feature, 1, function(x) mean(x != 0, na.rm = T))
    F.mtx = o.feature

    f.mtx[which(F.mtx == 0)] <- (1 - f.mtx[which(F.mtx == 0)])
    #index<-order(f.mtx)
    #f.mtx<-sort(f.mtx,decreasing = FALSE)
    F.mtx1<-F.mtx
    #for(i in length(index)){
    #  F.mtx1[i]<-F.mtx[which(index==i)]
    #}
    #num_subject<-round(nr*trim)
    #F.mtx1<-F.mtx1[num_subject:nr]
    #f.mtx<-f.mtx[num_subject:nr]
    feature_sentivity<-sum(f.mtx[which(F.mtx1 == 1)])/sum(F.mtx1 == 1)
    if(sum(F.mtx1==0)!=0){
      feature_specficity<-sum(f.mtx[which(F.mtx1 == 0)])/sum(F.mtx1 == 0)
      feature_score<-feature_sentivity+feature_specficity-1
    }else{
      feature_score<-feature_sentivity
    }



    #U.min2 = 2 - (mean(r.mtx[which(o.mtx == 1)],na.rm = T) + mean(r.mtx[which(o.mtx == 0)],na.rm =T))
    #U.min3 = 1 - mean(r.mtx[which(o.mtx == 1)],na.rm = T)
    final_score<-(clus_score+feature_score)/2
    return(list(final_score=final_score,clus_score=clus_score))
  }
  result_clus<-list(1)
  result_final<-list(1)
  for(i in 1:length(k_vector)){
    #set.seed(12315)
    temp<-lapply(wbounds_list[[i]], function(wbounds,x,k){return(score_wbounds(wbounds,x,k))},x,k_vector[i])
    result_clus[[i]]<-rep(-1,length(wbounds_list[[i]]))
    result_final[[i]]<-rep(-1,length(wbounds_list[[i]]))
    for(i1 in 1:length(wbounds_list[[i]])){
      result_clus[[i]][i1]<-temp[i1][[1]]$clus_score
      result_final[[i]][i1]<-temp[i1][[1]]$final_score
    }

  }

  for(i1 in 1:length(k_vector)){
    if(sum(is.na(result_clus[[i1]]))!=0){
      index_na<-which(is.na(result_clus[[i1]]))
      for(i2 in 1:length(index_na)){
        result_clus[[i1]][index_na[i2]]<-score_wbounds(wbounds_list[[i1]][index_na[i2]],x,k_vector[i1])$clus_score
        result_final[[i1]][index_na[i2]]<-score_wbounds(wbounds_list[[i1]][index_na[i2]],x,k_vector[i1])$final_score
      }
    }
  }
  result_k1<-rep(-1,length(k_vector))
  for(i in 1:length(k_vector)){
    result_k1[i]<-max(result_clus[[i]])
  }
  res_k1<-k_vector[max(which(result_k1==max(result_k1)))]
  index<-which.max(result_final[[which(res_k1==k_vector)]])
  res_w1<-wbounds_list[[which(res_k1==k_vector)]][index]



  matrix_length<-max(unlist(lapply(wbounds_list,length)))
  res_matrix<-matrix(rep(-1,length(k_vector)*matrix_length),nrow=length(k_vector))
  for(i in 1:length(k_vector)){
    length_temp<-length(result_final[[i]])
    res_matrix[i,1:length_temp]<-result_final[[i]]
  }


  res2<-which(res_matrix==max(res_matrix),arr.ind=T)
  index_temp<-which.max(res2[,1])
  res2<-res2[index_temp,]
  res_k2<-res2[1]
  res_w2<-wbounds_list[[res_k2]][res2[2]]
  res_k2<-k_vector[res_k2]
  #K.estimate = wbounds[which.max(final_score)]

  #plot(wbounds, final_score)

  return(list(res_k1=res_k1,res_k2=res_k2,res_w1=res_w1,res_w2=res_w2,result_clus=result_clus,result_final=result_final))
}
