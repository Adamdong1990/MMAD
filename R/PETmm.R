#' Poisson model for transmission tomography
#'
#' @param theta attenuation coefficient
#' @param y transmission measurement of detector
#' @param A the n*q system matrix
#' @param s the blank scan counts
#' @param r the mean number of background counts of detector
#'
#' @return requires a value
#' @export
#'
#' @examples
#' y=c(4,1); theta=c(3,-4);s=c(3.7, 3.8); r=c(0,0)
#' aa = c(0.029, 0.027, 0.007, 0.06)
#' A=matrix(aa,2)
#' PETmm(y, A, s, r, theta)
PETmm <- function(y, A, s, r, theta, Maxiter = 1000, convergence = 1e-06, method = "ADMM", ...)
{
  n = length(y)
  p = length(theta)
  TH = matrix(rep(theta,each=n),n,p)

  # log-likelihood function
  a = exp(-rowSums(A*TH))
  b = log(r+s*a)
  log_ell = sum(-s*a+y*b)
  el = c(log_ell)

  error = 3
  result <- list()

  for (k in 1:Maxiter)
  {
    if (error > convergence)
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
      print_err <- error
      print_k <- k
    }
  }

  Fisher_Matrix <- -n*PETFx(A, s, r, theta)
  std_err <- sqrt(diag(solve(Fisher_Matrix)))
  std_theta <- std_err

  # confidence intervals
  ci_theta_lower <- theta - 1.96*std_theta
  ci_theta_upper <- theta + 1.96*std_theta

  ELL = el[length(el)]
  Rate = PET_Rate(A, s, r, theta)

  #add values of AIC and BIC
  aic <- (2 * length(theta)) - (2 * ELL)
  bic <- log(length(y)) * length(theta) - 2 * ELL
  info_criteria <- c(AIC=aic, BIC=bic)

  result$call <- match.call()
  result$print_n <- n
  result$print_k <- print_k
  result$print_err <- print_err
  result$ELL <- ELL
  result$theta <- theta
  result$std_theta <- std_theta
  result$ci_theta_lower <- ci_theta_lower
  result$ci_theta_upper <- ci_theta_upper

  result$Rate <- Rate
  result$info_criteria <- info_criteria
  result$convergence <- convergence
  return(result)
}
