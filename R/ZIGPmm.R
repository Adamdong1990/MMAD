
#' MM algorithm based on the AD method for Multivariate compound zero-inflated generalized Poisson distribution
#'
#' @param phi0 probability
#' @param phi probability
#' @param la scale parameter
#' @param th discrete parameter
#' @param y matrix data
#'
#' @return
#' @export
#'
#' @examples
#' y = matrix(7,1,1);
#' ZIGPmm(0.5,0.5,1,0.1,y)
ZIGPmm <- function(y, phi0, phi, la, th, Maxiter = 500, convergence = 1e-06, method = "ADMM", ...)
{
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

  phi0_std <- c(phi0)
  phi_std <- c(phi)
  la_std <- c(la)
  th_std <- c(th)
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

      phi0_std <- append(phi0_std, phi0)
      phi_std <- rbind(phi_std, phi)
      la_std <- rbind(la_std, la)
      th_std <- rbind(th_std, th)

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

  std_phi0 <- sd(phi0_std)/sqrt(length(phi0_std))
  std_phi <- apply(phi_std, 2, sd)/sqrt(length(phi_std)/m)
  std_la <- apply(la_std, 2, sd)/sqrt(length(la_std)/m)
  std_th <- apply(th_std, 2, sd)/sqrt(length(th_std)/m)

  phi0_t_val <- phi0/std_phi0
  phi_t_val <- phi/std_phi
  la_t_val <- la/std_la
  th_t_val <- th/std_th

  ELL = el[length(el)]
  alpha = c(phi0, phi, la, th)
  Rate = ZIGP_CRate(phi0, phi, la, th)

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
  result$phi0_t_val <- phi0_t_val
  result$phi <- phi
  result$std_phi <- std_phi
  result$phi_t_val <- phi_t_val
  result$la <- la
  result$std_la <- std_la
  result$la_t_val <- la_t_val
  result$th <- th
  result$std_th <- std_th
  result$th_t_val <- th_t_val

  result$Rate <- Rate
  result$info_criteria <- info_criteria
  result$convergence <- convergence
  return(result)
}
