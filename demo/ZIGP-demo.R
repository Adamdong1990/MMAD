
rm(list=ls())
setwd("C:/code")


ZIGP_sample <- function(phi,la,th,n)
{
  GP <- function(x, la0, th0)
  {
    b<-gamma(x+1)
    if(th0 <1 && th0>=0) {
      p<-exp(-la0-th0*x)*la0*(la0+th0*x)^(x-1)/b }
    if(th0<0 && th0 >-1 &&  la0+th0*x <=0){
      p<-0 }
    if(th0<0 && th0>-1 && la0+th0*x > 0) {
      p<-exp(-la0-th0*x)*la0*(la0+th0*x)^(x-1)/b }
    return(p)
  }

  m <- length(la)
  X <- matrix(0,n,m)
  P <- matrix(0,100,m)
  for(j in 1:m)
  {
    for(i in 0:99)
    {   P[i+1,j] <- GP(i,la[j],th[j])     }
    a<-0:99
    X[,j]<- sample(a,n,P[,j],replace=TRUE)
  }
  z <- sample(c(0,1),n,c(phi,1-phi),replace=TRUE)
  zz <- matrix(z,n,m)
  y <- zz*X
  return(y)
}


N=500
RES = matrix(0,N,204)
for(t in 1:N)
{
  n = 300
  phi0 = 0.2
  la0 = rep(7,100)
  th0 = rep(0.3,100)
  y = ZIGP_sample(phi0,la0,th0,n)

  phi = 0.5
  la = rep(1,100)
  th = rep(0.1,100)

  result <- ZIGPmm(phi, la, th, y)
  RES[t,] = c(result$k, result$ELL, result$alpha, result$mse)


}

MRES = apply(RES,2,mean)
MRES
