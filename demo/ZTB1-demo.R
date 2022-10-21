
rm(list=ls())
setwd("C:/code")

ZTB_sample <- function(m,th,n)
{
  x <- rep(0,n)
  for(i in 1:n)
  {
    repeat
    {  a <- rbinom(1,m,th)
    if(a>0)   break}
    x[i] <- a
  }
  return(x)
}


N=500
RES = matrix(0,N,4)
for(i in 1:N)
{
  n = 200
  m = 10
  th0 = 0.6
  y = ZTB_sample(m,th0,n)

  th = 0.1
  result = ZTB1mm(m,th,y)

  RES[i,] = c(result$k, result$ELL, result$th, result$mse)

}
MRES = apply(RES,2,mean)
MRES
