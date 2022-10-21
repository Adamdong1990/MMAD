#' MM algorithm based on the AD method for Multivariate compound zero-inflated generalized Poisson distribution
#'
#' @param phi probability
#' @param la scale parameter
#' @param th discrete parameter
#' @param y matrix data
#'
#' @return data
#' @export
#'
#' @examples
#' y = matrix(7,1,1); ZIGPmm(0.5,1,0.1,y)
ZIGPmm <- function(phi, la, th, y)
{
  n = nrow(y)
  m = length(la)
  zero = matrix(0,n,m)
  yz = apply(1*(y==zero),1,prod)
  n0 = sum(yz)

  # log-likelihood function
  a = n0*log(phi+(1-phi)*exp(-sum(la))) + (n-n0)*log(1-phi) - (n-n0)*sum(la)
  b = rep(0,m)
  for(i in 1:m)
  {
    b[i] = n*log(la[i])-sum(y[,i])*th[i] + sum((y[,i]-1)*log(la[i]+th[i]*y[,i]))
  }
  log_ell = a + sum(b)
  el = c(log_ell)

  k = 1
  error = 3
  alpha0 = c(phi,la,th)

  while( error > 1e-06 )
  {
    beta = phi + (1-phi)*exp(-sum(la))
    phi = n0*phi/(n*beta)

    for(i in 1:m)
    {
      la[i] = ( n+sum( (y[,i]-1)*la[i]/(la[i]+ th[i]*y[,i]) ) )/(n-n*phi)
      th[i] = sum( th[i]*y[,i]*(y[,i]-1)/(la[i]+th[i]*y[,i]) )/sum(y[,i])

      b[i] = n*log(la[i])-sum(y[,i])*th[i] + sum((y[,i]-1)*log(la[i]+th[i]*y[,i]))
    }

    a = n0*log(phi+(1-phi)*exp(-sum(la))) + (n-n0)*log(1-phi) - (n-n0)*sum(la)

    log_el = a + sum(b)
    el <- append(el, log_el)
    error = abs(el[k+1]-el[k])/(1+abs(el[k]))
    k = k + 1
  }

  ELL = el[length(el)]
  alpha = c(phi,la,th)
  mse = sum( (alpha-alpha0)^2 )/(2*m+1)

  result <- list()
  result$k <- k
  result$ELL <- ELL
  result$alpha <- alpha
  result$mse <- mse
  return(result)
}
