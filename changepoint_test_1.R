# setwd("C:/Users/jd033/Google Drive/Change Points/R")

rm(list = ls())
setwd("C:/Users/jd033/OneDrive/Documents/R/changepoint")

cut<-sort(runif(3))

cut<-c(.3,.7)
n<-250
Y<-rep(0,n)
for(i in 1:n){
  if(i/n <= .3){Y[i]<-rnorm(1)}
  if(i/n > .3  && i/n <=.7){Y[i]<-rnorm(1)+2.5}
  
  if(i/n > .7){Y[i]<-rnorm(1)-2.5}
  
}
plot(Y)

source("est_cp.R")

fit<- est_cp(Y)

N = dim(fit$estmat)[2]
  
yhat <-  apply(fit$estmat[,(N-1000):N],1,mean)
points(yhat,type="l")

chp<-apply(fit$indmat[,4000:5000],1,mean)
plot(chp)
