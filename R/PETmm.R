#' Poisson model for transmission tomography
#'
#' @param theta attenuation coefficient
#' @param y transmission measurement of detector
#' @param A the n*q system matrix
#' @param s the blank scan counts
#' @param r the mean number of background counts of detector
#'
#' @return
#' @export
#'
#' @examples
#' y=c(4,1); theta=c(3,-4);s=c(3.7, 3.8); r=c(0,0)
#' aa = c(0.029, 0.027, 0.007, 0.06)
#' A=matrix(aa,2)
#' PETmm(theta,y,A,s,r)
PETmm <- function(theta,y,A,s,r)
{
  k=1
  error = 3
  theta0 = theta
  n = length(y)
  p = length(theta)
  TH = matrix(rep(theta,each=n),n,p)
  a = exp(-rowSums(A*TH))
  b = log(r+s*a)

  # log-likelihood function
  log_ell = sum(-s*a+y*b)
  el = c(log_ell)

  while( error>1e-06 )
  {
    ww = t(matrix(rep(rowSums(A),each=p),p,n))
    w = A/ww
    for(j in 1:p)
    {
      Q1 = sum(A[,j]*s*a)-sum(y*A[,j]*s*a/(r+s*a))
      Q2 = -sum(A[,j]^2*s*a/w[,j])
      theta[j] = theta[j]-Q1/Q2
    }

    TH = matrix(rep(theta,each=n),n,p)
    a = exp(-rowSums(A*TH))
    b = log(r+s*a)
    log_el = sum(-s*a+y*b)
    el <- append(el, log_el)
    error = abs(el[k+1]-el[k])/(1+abs(el[k]))
    k=k+1
  }

  ELL = el[length(el)]
  mse = sum( (theta-theta0)^2 )/(2*p+1)

  result <- list()
  result$k <- k
  result$ELL <- ELL
  result$theta <- theta
  result$mse <- mse
  return(result)
}
