rm(list=ls())


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


a = 5; n = 1000; mu = 7; si2 = 4


  y = LTN_sample(a, mu, si2, n)

  result <- LTNmm(y, 5, 1, 1)
  summary.LTN(result)


  devtools::install_github("williazo/tcensReg")
  library(tcensReg)

  dt <- data.frame(y)
  t_mod <- tcensReg(y ~ 1, data=dt, a=5)
  summary(t_mod)

  mu = 7.1197; si2 = 3.7988
  a1 = (a-mu)/(sqrt(si2))
  dphi = -a1*exp(-0.5*a1^2)/sqrt(2*pi)
  c = 1-pnorm(a1)

  -1/(2*si2^2) + (a1^2*dphi-a1*dnorm(a1))/(4*c*si2^2)+ (4+a1^2)*(dnorm(a1))^2/(4*c^2*si2^2)


