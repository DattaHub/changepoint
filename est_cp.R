est_cp <- function(Y, niter = 5000, verbose = TRUE){
  
  n = length(Y)
  indvec<-c(1:n)-c(1:n)
  pos<-c(0)
  indmat<-matrix(0,n,niter)
  estmat<-indmat
  
  pr<-1/n
  
  for(it in 1:niter)
  {
    if(isTRUE(verbose) && it %% 10 == 0)
      cat("Iteration ",it, "\n")
    
    for(i in sample(2:(n-1)))
      
    {
      tempind<-indvec
      posold<-pos
      
      l<-max(pos[which(pos<i)])
      m<-min(pos[which(pos>i)],(n+1))
      
      Y1<-Y[(l+1):(i-1)]
      Y2<-Y[(i):(m-1)]
      
      Y3<-Y[(l+1):(m-1)]
      
      
      ss1<-sum((Y1-mean(Y1))^2)/2
      ss2<-sum((Y2-mean(Y2))^2)/2
      ss3<-sum((Y3-mean(Y3))^2)/2
      
      l1<-(i-l)
      l2<-m-i
      
      delta<-ss3-ss2-ss1+0.5*(log(l1)+log(l2)-log(l1+l2))
      
      ratio<-delta+log(pr)-log(1-pr)
      eratio<-exp(ratio)/(1+exp(ratio))
      u<-runif(1)
      
      tempind[i]<-0
      if(u<eratio){ tempind[i]<-1}
      
      indvec<-tempind
      posnew<-c(0,c(which(indvec==1)))
      
      pos<-posnew
 
      # l<-max(pos[which(pos<i)])
      # m<-min(pos[which(pos>i)],n)

      estmat[i,it]<-mean(Y2)
      if(indvec[i]==0){estmat[i,it]<-mean(Y3)}

    }
    estmat[n,it]<-mean(Y[max(pos):n])
    indmat[,it]<-indvec
  }
return(list(estmat = estmat, indmat = indmat))
}