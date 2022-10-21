rm(list=ls())
setwd("C:/code")


LTN_sample <- function(a, mu, si2, n)
{
  x <- rep(0,n)
  for (i in 1:n)
  {
    repeat {
      z <- rnorm(1,mu,sqrt(si2))
      if(z >= a)
      {break} }
    x[i] <- z
  }
  return(x)
}


N=500; a = 5; n = 200; mu = 7; si2 = 4

RES = matrix(0,N,5)
for(i in 1:N)
{
  y = LTN_sample(a, mu, si2, n)
  result <- LTNmm(a, mu, si2, y)
  RES[i,] = c(result$k, result$ELL, result$alpha, result$mse)

}

MRES = apply(RES,2,mean)
MRES

