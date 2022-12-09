
#' MM algorithm based on the AD method for Multivariate compound zero-inflated generalized Poisson distribution
#'
#' @param phi0 probability
#' @param phi probability
#' @param la scale parameter
#' @param th discrete parameter
#' @param y matrix data
#'
#' @return requires a value
#' @export
#'
#' @examples
#' y = matrix(7,1,1);
#' CZIGPmm(0.5,0.5,1,0.1,y)
CZIGPmm <- function(y, phi0, phi, la, th, Maxiter = 500, convergence = 1e-06, method = "ADMM", ...)
{
  if(any(phi0 < 0) | any(phi0 > 1))
  { stop("The probability `phi0` of Bernoulli distribution should be between 0 and 1", call.=FALSE) }

  if(any(phi < 0) | any(phi > 1))
  { stop("The probability `phi` of Bernoulli distribution should be between 0 and 1", call.=FALSE) }

  if(length(phi) == 0 | length(la) == 0 | length(th) == 0)
  { stop("The length of the variable `phi` or `la` or `th` is not allowed to be 0", call.=FALSE) }

  if(length(phi) != length(la))
  { stop("The length of the variable `phi` and `la` must be the same", call.=FALSE) }
  if(length(la) != length(th))
  { stop("The length of the variable `la` and `th` must be the same", call.=FALSE) }

  n = nrow(y)
  m = length(la)
  zero = matrix(0,n,m)
  yz = apply(1*(y==zero),1,prod)
  n0 = sum(yz)
  Iy = 1-yz

  # log-likelihood function
  a = n0*log( phi0+(1-phi0)*prod(phi+(1-phi)*exp(-la)) ) + (n-n0)*log(1-phi0)
  b1 = sum(Iy*(y==0)*log(phi+(1-phi)*exp(-la)))
  b2 = sum( Iy*(y!=0)*(log(1-phi)+log(la)+(y-1)*log(la+th*y)-la-th*y) )
  log_ell = a + b1 + b2
  el = c(log_ell)

  error = 3
  alpha0 = c(phi0,phi,la,th)
  result <- list()

  for (k in 1:Maxiter)
  {
    if (error > convergence)
    {
      be0 = phi0+(1-phi0)*prod(phi+(1-phi)*exp(-la))
      phi0 = n0*phi0/(n*be0)
      be = phi+(1-phi)*exp(-la)
      for(i in 1:m)
      {
        ap = n0*phi[i]*(be0-phi0)/(be0*be[i]) + sum(Iy*(y[,i]==0)*phi[i]/be[i])
        phi[i] = ap/(n-n0*phi0/be0)
        ala = sum(Iy*(y[,i]!=0)*(la[i]+th[i])*y[,i]/(la[i]+th[i]*y[,i]))
        bla = n0*(be0-phi0)*(be[i]-phi[i])/(be0*be[i]) + sum( Iy*(1-(y[,i]==0)*phi[i]/be[i]) )
        la[i] = ala/bla
        at = sum( Iy*(y[,i]!=0)*th[i]*y[,i]*(y[,i]-1)/(la[i]+th[i]*y[,i]) )
        bt = sum(Iy*(y[,i]!=0)*y[,i])
        th[i] = at/bt
      }

      a = n0*log( phi0+(1-phi0)*prod(phi+(1-phi)*exp(-la)) ) + (n-n0)*log(1-phi0)
      b1 = sum(Iy*(y==0)*log(phi+(1-phi)*exp(-la)))
      b2 = sum( Iy*(y!=0)*(log(1-phi)+log(la)+(y-1)*log(la+th*y)-la-th*y) )
      log_el = a + b1 + b2
      el <- append(el, log_el)
      error=abs(el[k+1]-el[k])/(1+abs(el[k]))

      print_err <- error
      print_k <- k
    }
  }

  Fisher_Matrix <- -n*CZIGPFx(phi0, phi, la, th)
  std_err <- sqrt(diag(solve(Fisher_Matrix)))
  std_phi0 <- std_err[1]
  std_phi <- std_err[2: (1+m)]
  std_la <- std_err[(2+m): (1+2*m)]
  std_th <- std_err[(2+2*m): (1+3*m)]

  # confidence intervals
  ci_phi0_lower <- phi0 - 1.96*std_phi0
  ci_phi0_upper <- phi0 + 1.96*std_phi0

  ci_phi_lower <- phi - 1.96*std_phi
  ci_phi_upper <- phi + 1.96*std_phi

  ci_la_lower <- la - 1.96*std_la
  ci_la_upper <- la + 1.96*std_la

  ci_th_lower <- th - 1.96*std_th
  ci_th_upper <- th + 1.96*std_th


  ELL = el[length(el)]
  alpha = c(phi0, phi, la, th)
  Rate = CZIGP_CRate(phi0, phi, la, th)

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

  result$phi <- phi
  result$std_phi <- std_phi
  result$ci_phi_lower <- ci_phi_lower
  result$ci_phi_upper <- ci_phi_upper

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
