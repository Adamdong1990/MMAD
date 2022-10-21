rm(list=ls())
setwd("C:/code")



LG_sample = function(n,theta)
{
  p = length(theta)
  A = matrix(0,n,p)
  r = rpois(n,0.1)
  z = rnorm(n,0,1)
  s = exp(1+0.3*z)
  y = rep(0,n)
  for(i in 1:n)
  {
    A[i,] = runif(p,min=0,max=0.1)
    la = r[i] + s[i]*exp(-sum(A[i,]*theta))
    y[i] = rpois(1, la)
  }
  return(list(y=y,A=A,s=s,r=r))
}


N = 500
RES = matrix(0,N,33)
for(i in 1:N)
{
  n = 500
  theta1=rep(0.5,10)
  theta2=rep(0.8,10)
  theta3=rep(0.2,10)
  theta = c(theta1,theta2,theta3)
  p = length(theta)
  yy = LG_sample(n,theta)

  y = yy$y; A = yy$A; s = yy$s; r = yy$r
  theta = rep(0.1,p)

  result <- PETmm(theta,y,A,s,r)
  RES[i,] = c(result$k, result$ELL, result$theta, result$mse)

}

MRES = apply(RES,2,mean)
MRES
