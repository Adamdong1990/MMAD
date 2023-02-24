# PETFx
PETFx <- function(A, s, r, theta)
{
  n = length(s)
  p = length(theta)
  df = matrix(0, p, p)

  TH = matrix(rep(theta,each=n),n,p)
  a = exp(-rowSums(A*TH))

  for(i in 1:p)
  {
    for(j in 1:p)
    {
      df[j,i] = df[i,j] = sum( -A[,j]*A[,i]*s*a + s*r*a*A[,j]*A[,i]/(r+s*a) )
    }
  }

  return(df)
}

# PETQx
PETQx <- function(A, s, r, theta)
{
  n = length(s)
  p = length(theta)
  dth = c()

  TH = matrix(rep(theta,each=n),n,p)
  a = exp(-rowSums(A*TH))
  ww = t(matrix(rep(rowSums(A),each=p),p,n))
  w = A/ww

  for(j in 1:p)
  {
    th <- -sum(A[,j]^2*s*a/w[,j])
    dth <- append(dth, th)
  }

  dQ = diag(dth)
  return(dQ)
}

# PET_Rate
PET_Rate <- function(A, s, r, theta)
{
  DF = PETFx(A, s, r, theta)
  DGI = solve(PETQx(A, s, r, theta))
  A =  DGI%*%DF
  rate = 1-min(eigen(A)$values)
  return(rate)
}
