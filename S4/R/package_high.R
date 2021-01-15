
# Output = vector("list", 10)
# for(i in 1:10){
#   error = F
#   while(!error){
#     output = tryCatch(functionOfLambda(lambda[i]), error = function(x) NA)
#     if(is.na(output)) error = T
#   }
#   Output[[i]] = output
# }


########################put parallel option: prediction strength


score_PS<-function(x,wbounds,K,cv){
  n = nrow(x)
  x_fold <- createFolds(1:n, k = cv, list = TRUE, returnTrain = FALSE)
  result.subsampling = lapply(1:cv, function(b){
    #x_fold <- createFolds(1:n, k = cv, list = TRUE, returnTrain = FALSE)

    index_test<-x_fold[[b]]
    index_train<-unlist(x_fold[-b])
    train_x<-x[as.numeric(index_train),]
    test_x<-x[index_test,]

    #km.out <- kmeans(xb, centers = K, nstart = 100)
    #hc = hclust(dist(xb), method = "complete", members = NULL)
    #cl<-cmeans(xb,3,m=K)
    te = KMeansSparseCluster1(test_x, K=K, wbounds = wbounds, nstart = 20, ##########
                              silent = T, maxiter=6, centers=NULL)
    tr=KMeansSparseCluster1(train_x, K=K, wbounds = wbounds, nstart = 20, ##########
                            silent = T, maxiter=6, centers=NULL)
    select_feature_train<-as.numeric(te[[1]]$ws>0)
    select_feature_test<-as.numeric(tr[[1]]$ws>0)
    cluster_train<-tr[[1]]$Cs
    cluster_test<-te[[1]]$Cs
    index_choose<-which(select_feature_test==1)
    final_feature<-sum(select_feature_train[index_choose])/length(index_choose)
    centroid<-matrix(rep(0,K*dim(x)[2]),nrow=K)
    ######calculate centroid
    for(i in 1:K){
      if(sum(cluster_train==i)==1){
        centroid[i,]<-train_x[which(cluster_train==i),]
        next
      }
      centroid[i,]<-apply(train_x[which(cluster_train==i),],2,mean)
    }
    ############assign test data class label according to train data
    predict_test<-rep(0,nrow(test_x))
    for(i in 1:length(predict_test)){
      temp<-rep(0,K)
      for(i1 in 1:K){
        temp[i1]<-weighted_l2(centroid[i1,],test_x[i,],tr[[1]]$ws)
      }
      predict_test[i]<-which.min(temp)
    }
    # predict_test<-rep(0,nrow(test_x))
    # for(i in 1:length(predict_test)){
    #   temp<-rep(0,k)
    #   for(i1 in 1:k){
    #     temp[i1]<-dist(rbind(centroid[i1],test_x[i]))
    #   }
    #   predict_test[i]<-which.max(temp)
    # }
    ##############calculate prediction score for each cluster
    score_eachcluster<-rep(-1,K)
    for(i in 1:K){
      n1<-length(which(cluster_test==i))
      temp1<-predict_test[which(cluster_test==i)]
      consensus.matrix = sapply(1:length(temp1), function(i){
        as.integer(temp1[i] == temp1)
      })
      if(n1==1){
        score_eachcluster[i]<-1
        next
      }
      diag(consensus.matrix)<-0
      score_eachcluster[i]<-sum(consensus.matrix)/(n1*(n1-1))
    }
    final_clus<-min(score_eachcluster)
    ##################the prediction for the feature
    #sel_feature_test<-which(te[[1]]$ws>0)
    #final_feature<-length(intersect(sel_feature_train,sel_feature_test))/length(sel_feature_test)

    #final_score<-(final_clus+final_feature)/2
    final_score<-(final_clus+final_feature)/2

    return(list(final_score = final_score,clus_score=final_clus))
  })

  res.final = sapply(result.subsampling, function(x) x$final_score, simplify = "array")
  res.clus = sapply(result.subsampling, function(x) x$clus_score, simplify = "array")
  return(list(final=mean(res.final),clus=mean(res.clus)))
}


######put parallel option S4

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
  stat<-sort(stat,decreasing = T)
  #plot(stat)
  #title(paste("plot of S4.Iter","_K=",K,sep=""))
  #barplot(stat,width = 0.832)
  #title(paste("k=",K,sep=""))
  #abline(v=n*0.95,col="red",lwd=3.5)
  clus_score = mean(stat)
  #clus_score = mean(stat[1:round(nr*per)])


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


# S4_highdim = function(x, lambda_list,trim=0.05,k_vector=c(2,3,4),n.resample=25,num.cores=1){
#   result_clus<-list(1)
#   result_final<-list(1)
#   for(ind_k in 1:length(k_vector)){
#     k<-k_vector[ind_k]
#     lambda_vector<-lambda_list[[ind_k]]
#     wcs = mclapply(1:length(lambda_vector),function(ind_lambda,k,x){
#       lambda<-lambda_vector[ind_lambda]
#       error = T
#       ind<-0
#       while(error){
#         result= tryCatch(score_lambda(lambda,x,k,trim=0.05,n.resample=25), error = function(x) NA)
#         if(!is.na(result[1])){
#           error = F
#         }
#         ind<-ind+1
#         if(ind>3){break}
#       }
#
#
#       #result<-score_lambda(lambda,x,k,trim=0.05,n.resample=25)
#       return(result)
#     },k = k,x=x,mc.cores = num.cores)
#     result_clus[[ind_k]]<-NA
#     result_final[[ind_k]]<-NA
#     for(ind_lambda in 1:length(lambda_vector)){
#       result_clus[[ind_k]][ind_lambda]<-wcs[[ind_lambda]]$clus_score
#       result_final[[ind_k]][ind_lambda]<-wcs[[ind_lambda]]$final_score
#     }
#     }
#   result_k<-rep(-1,length(k_vector))
#   for(i in 1:length(k_vector)){
#     result_k[i]<-max(result_clus[[i]])
#   }
#   res_k<-k_vector[max(which(result_k==max(result_k)))]
#   index<-which.max(result_final[[which(res_k==k_vector)]])
#   res_lambda<-lambda_list[[which(res_k==k_vector)]][index]
#   return(list(clus_score=result_clus,final_score=result_final,optimal_k=res_k,optimal_lambda=res_lambda))
# }

# j<-1
# k<-1
# q = c(50,200) #number of DE features out of 1000
# u = c(0.8,0.6,0.4)
# h<-100#DE evidence
# set.seed(i)
# x<-simulation_data_simple(h=h,q=q[j],u=u[k])
# lambda_list<-list(1)
# k_vector<-c(2,3)
# for(l in 1:length(k_vector)){
#   lambda_list[[l]] = region.lamda(b=10,x=x,h=h,K=k_vector[l])
# }
# S4_highdim(x, wbounds_list,trim=0.05,k_vector=c(2,3),n.resample=25,num.cores=2)
#
#################################gap statistic function
# gap_one<-function(x,K=3,wbounds,n.perms){
#   a = KmeansSparseCluster.permute1(x, K=K, nperms = n.perms, wbounds = wbounds,
#                                    silent = FALSE, nvals = 20, centers=NULL)
#   maxgap_minus_sd<-max(a$gaps)-a$sdgaps[which.max(a$gaps)]
#   return(wbounds[which.max(a$gaps>maxgap_minus_sd)])
# }

# gap_two<-function(x,k_vector=c(2,3,4),wbounds_list,n.perms=25,num_cores=1){
#   res_gap<-list(1)
#   res_gap_sd<-list(1)
#   res = mclapply(1:length(k_vector),function(l){
#     a = KmeansSparseCluster.permute1(x, K=k_vector[l], nperms = n.perms, wbounds = wbounds_list[[l]],
#                                      silent = FALSE, nvals = 20, centers=NULL)
#     return(list(gap=a$gaps,sd=a$sdgaps))
#   },mc.cores = num_cores)
#
#
#   for(i in 1:length(res)){
#     res_gap[[i]]<-res[[i]]$gap
#     res_gap_sd[[i]]<-res[[i]]$sd
#   }
#
#
#   maxgap_k<-rep(-1,length(k_vector))
#   for(i in 1:length(k_vector)){
#     maxgap_k[i]<-max(res_gap[[i]])
#   }
#   gap_k<-k_vector[max(which(max(maxgap_k)==maxgap_k))]
#   #a<-KmeansSparseCluster.permute1(x, K=gap_k, nperms = n.perms, wbounds = wbounds_list[[which(k_vector==gap_k)]],
#   #                                silent = FALSE, nvals = 20, centers=NULL)
#   #maxgap_minus_sd<-max(a$gaps)-a$sdgaps[which.max(a$gaps)]
#   a<-res_gap[[which(k_vector==gap_k)]]
#   maxgap_minus_sd<-max(a)-res_gap_sd[[which(k_vector==gap_k)]][which.max(a)]
#   res_w<-wbounds_list[[which(k_vector==gap_k)]][which.max(a>maxgap_minus_sd)]
#   return(list(res_k=gap_k,res_w=res_w))
# }

#############################simulation of data
standardize_cov<-function(x){
  y<-x
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(x)[1]){
      y[i,j]<-x[i,j]/(sqrt(x[i,i])*sqrt(x[j,j]))
    }
  }
  return(y)
}





simulation_data_cov<-function(K,POIN1=POIN1,POIN2=POIN2,POIN3=POIN3,Nm=POIM,M,Ua,Ub,Ulow,Uupper,sigma1,sigma2,Nnoise=600,cov=0.5){
  ##sigma1 is the biological variation of predictive genes
  ##sigma2 is the bilogical variation of the noise genes
  N1<-rpois(1,POIN1)
  N2<-rpois(1,POIN2)
  N3<-rpois(1,POIN3)
  N<-N1+N2+N3
  indicator_sample<-rep(1:3,c(N1,N2,N3))
  num_genes_in_module<-rep(0,M)
  for(i in 1:length(num_genes_in_module)){
    num_genes_in_module[i]<-rpois(1,Nm)
    while(num_genes_in_module[i]==0){
      num_genes_in_module[i]<-rpois(1,Nm)
    }
  }
  #num_genes_in_module<-rpois(M,Nm)
  indicator_module<-rep(1:M,num_genes_in_module)
  template_expression<-matrix(rep(0,M*K),nrow=K)
  for(i in 1:M){
    template_expression[,i]<-runif(K,Ua,Ub)
    while(((max(template_expression[,i])-min(template_expression[,i]))<Ulow)|
          ((max(template_expression[,i])-min(template_expression[,i]))>Uupper)){
      template_expression[,i]<-runif(K,Ua,Ub)
    }

  }
  #apply(template_expression,2,function(x){return(max(x)-min(x))})
  #####add biological variation
  gene_expression<-matrix(rep(0,N*M),nrow=N)
  for(i in 1:M){
    for(j in 1:N){
      gene_expression[j,i]<-rnorm(1,template_expression[indicator_sample[j],i],sigma1)
    }
  }
  #apply(gene_expression[which(indicator_sample==1),],2,mean)
  data<-matrix(rep(0,N*length(indicator_module)),nrow=N)
  for(i in 1:N){
    for(m in 1:M){
      v1<-riwish(60,((1-cov)*diag(num_genes_in_module[m])+
                       cov*matrix(rep(1,num_genes_in_module[m]^2),ncol=num_genes_in_module[m])))
      v<-standardize_cov(v1)
      data[i,which(indicator_module==m)]<-rmvnorm(1,mean=rep(gene_expression[i,m],num_genes_in_module[m]),sigma=v)
    }
  }
  #######add housekeeping genes
  gene_expression_template_housekeeping<-runif(Nnoise,Ua,Ub)
  for(i in 1:Nnoise){
    data<-cbind(data,rnorm(nrow(data),mean=gene_expression_template_housekeeping[i],sigma2))
  }

  true_label<-c(rep(1,N1),rep(2,N2),rep(3,N3))
  true_feature<-c(rep(1,sum(num_genes_in_module)),rep(0,Nnoise))
  return(list(data=data,true_label=true_label,true_feature=true_feature))

}

# simulation_data_simple<-function(h,q,u){
#   #q = c(50,200) #number of DE features out of 1000
#   #u = c(1,0.8,0.6,0.4)
#   #h<-1000#DE evidence
#
#   x = rbind(matrix(rnorm(33*h),ncol = h),
#             matrix(rnorm(33*h),ncol = h),
#             matrix(rnorm(33*h),ncol = h))
#   x[1:33,1:q] <- x[1:33,1:q]-u
#   x[34:66,1:q] <- x[34:66,1:q]
#   x[67:99,1:q] <- x[67:99,1:q]+u
#   return(x)
# }
###############################################
