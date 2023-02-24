
rm(list=ls())


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

n = 500
theta1=rep(0.5,3)
theta2=rep(0.8,3)
theta3=rep(0.2,3)
theta = c(theta1,theta2,theta3)
p = length(theta)
yy = LG_sample(n,theta)

y = yy$y; A = yy$A; s = yy$s; r = yy$r
theta = rep(0.1,p)

result <- PETmm(y, A, s, r, theta)
summary.PET(result)
