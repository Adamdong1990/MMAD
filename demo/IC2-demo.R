rm(list=ls())
setwd("C:/code")


sample <- function(la,n)
{
  u = runif(n,0,1)
  w = runif(n,0,1)
  U = -log(u)/la
  W = -log(w)/la
  L = pmin(U,W)
  R = pmax(U,W)
  return(list(L=L,R=R))
}


N=100
RES = matrix(0,N,1004)
for(j in 1:N)
{
  n = 500
  y=sample(0.5,n)
  L = y$L
  R = y$R

  result = IC2mm(L, R)
  RES[j,] = c(result$k, result$ELL, result$p, result$M, result$AM)
}

MRES = apply(RES,2,mean)
MRES
