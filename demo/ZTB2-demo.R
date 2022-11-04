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



  n = 200
  m = 10
  th0 = 0.6
  y = ZTB_sample(m,th0,n)

  th = 0.1
  result = ZTB2mm(y, m, th)
  summary.ZTB(result)

