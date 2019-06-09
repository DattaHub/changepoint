
cut<-sort(runif(3))

cut<-c(.3,.7)
n<-250
y<-rep(0,n)
for(i in 1:n)

{
if(i/n <= .3){y[i]<-rnorm(1)}
if(i/n > .3  && i/n <=.7){y[i]<-rnorm(1)+2.5}

if(i/n > .7){y[i]<-rnorm(1)-2.5}

}
Y<-y


N<-5000
indvec<-c(1:n)-c(1:n)
pos<-c(0)
indmat<-matrix(0,n,N)
estmat<-indmat

pr<-1/n

for(it in 1:N)
{
  
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
    
fit<-(apply(estmat[,4000:N],1,mean))


plot(Y)
points(fit,type="l")

chp<-apply(indmat[,4000:5000],1,mean)
plot(chp)






