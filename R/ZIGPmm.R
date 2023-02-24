#' Title
#'
#' @param y matrix data
#' @param phi0 probability
#' @param la scale parameter
#' @param th discrete parameter
#' @param Maxiter The maximum number of iterations
#' @param convergence Convergence value
#' @param method optimized method
#' @param ...
#'
#' @return requires a value
#' @export
#'
#' @examples
#' y = matrix(7,1,1);
#' ZIGPmm(y,0.5,1,0.1)
ZIGPmm <- function(y, phi0, la, th, Maxiter = 500, convergence = 1e-06, method = "ADMM", ...)
{
  if(any(phi0 < 0) | any(phi0 > 1))
  { stop("The probability `phi0` of Bernoulli distribution should be between 0 and 1", call.=FALSE) }

  if(length(phi0) == 0 | length(la) == 0 | length(th) == 0)
  { stop("The length of the variable `phi0` or `la` or `th` is not allowed to be 0", call.=FALSE) }

  if(length(la) != length(th))
  { stop("The length of the variable `la` and `th` must be the same", call.=FALSE) }

  n = nrow(y)
  m = length(la)
  zero = matrix(0,n,m)
  yz = apply(1*(y==zero),1,prod)
  n0 = sum(yz)
  Iy = 1-yz

  # log-likelihood function
  a = n0*log( phi0+(1-phi0)*exp(-sum(la))) + (n-n0)*log(1-phi0)
  b1 = sum(Iy*(y==0)*la)
  b2 = sum( Iy*(y!=0)*(log(la)+(y-1)*log(la+th*y)-la-th*y) )
  b3 = sum( Iy*(y!=0)*log(factorial(y[y!=0])) )
  log_ell = a + b1 + b2 - b3
  el = c(log_ell)

  error = 3
  result <- list()

  for (k in 1:Maxiter)
  {
    if (error > convergence)
    {
      be0 = phi0+(1-phi0)*exp(-sum(la))
      phi0 = n0*phi0/(n*be0)
      for(i in 1:m)
      {
        ala = sum(Iy*(y[,i]!=0))+sum( Iy*(y[,i]!=0)*la[i]*(y[,i]-1)/(la[i]+th[i]*y[,i]) )
        bla = n*(1-phi0)
        la[i] = ala/bla
        at = sum( Iy*(y[,i]!=0)*th[i]*y[,i]*(y[,i]-1)/(la[i]+th[i]*y[,i]) )
        bt = sum(Iy*(y[,i]!=0)*y[,i])
        th[i] = at/bt
      }

      a = n0*log( phi0+(1-phi0)*exp(-sum(la))) + (n-n0)*log(1-phi0)
      b1 = sum(Iy*(y==0)*la)
      b2 = sum( Iy*(y!=0)*(log(la)+(y-1)*log(la+th*y)-la-th*y) )
      b3 = sum( Iy*(y!=0)*log(factorial(y[y!=0])) )
      log_el = a + b1 + b2 - b3
      el <- append(el, log_el)
      error=abs(el[k+1]-el[k])/(1+abs(el[k]))

      print_err <- error
      print_k <- k
    }
  }

  Fisher_Matrix <- -n*ZIGPFx(phi0, la, th)
  std_err <- sqrt(diag(solve(Fisher_Matrix)))
  std_phi0 <- std_err[1]
  std_la <- std_err[2: (1+m)]
  std_th <- std_err[(2+m): (1+2*m)]

  # confidence intervals
  ci_phi0_lower <- phi0 - 1.96*std_phi0
  ci_phi0_upper <- phi0 + 1.96*std_phi0

  ci_la_lower <- la - 1.96*std_la
  ci_la_upper <- la + 1.96*std_la

  ci_th_lower <- th - 1.96*std_th
  ci_th_upper <- th + 1.96*std_th


  ELL = el[length(el)]
  alpha = c(phi0, la, th)
  Rate = ZIGP_Rate(phi0, la, th)

  #add values of AIC and BIC
  aic <- (2 * length(alpha)) - (2 * ELL)
  bic <- log(length(y)) * length(alpha) - 2 * ELL
  info_criteria <- c(AIC=aic, BIC=bic)

  result$call <- match.call()
  result$print_n <- n
  result$print_k <- print_k
  result$print_err <- print_err
  result$ELL <- ELL
  result$phi0 <- phi0
  result$std_phi0 <- std_phi0
  result$ci_phi0_lower <- ci_phi0_lower
  result$ci_phi0_upper <- ci_phi0_upper

  result$la <- la
  result$std_la <- std_la
  result$ci_la_lower <- ci_la_lower
  result$ci_la_upper <- ci_la_upper

  result$th <- th
  result$std_th <- std_th
  result$ci_th_lower <- ci_th_lower
  result$ci_th_upper <- ci_th_upper

  result$Rate <- Rate
  result$info_criteria <- info_criteria
  result$convergence <- convergence
  return(result)
}
