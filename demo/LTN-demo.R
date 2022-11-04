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


N=500; a = 5; n = 100; mu = 7; si2 = 4

RES = matrix(0,N,5)
for(i in 1:N)
{
  y = LTN_sample(a, mu, si2, n)
  mu = 1; si2 = 1
  result <- LTNmm(y, 5, 1, 1)
  summary.LTNmm(result)
  RES[i,] = c(result$k, result$ELL, result$alpha, result$mse)

}

MRES = apply(RES,2,mean)
MRES

