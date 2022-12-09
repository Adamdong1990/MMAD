#' MM algorithm based on the AD method for Zero-truncated binomial distribution(1)
#'
#' @param m number of experiments
#' @param th probability of success
#' @param y requires a value
#'
#' @return requires a value
#' @export
#'
#' @examples
#' y = c(1, 2)
#' ZTB1mm(5, 0.2, y)
ZTB1mm <- function(y, m, th, Maxiter = 500, convergence = 1e-06, method = "ADMM", ...)
{
  # check 'm'
  if (!is.numeric(m))
  { stop("`m` must be numeric", call. = FALSE) }
  if (length(m) !=1)
  { stop("`m` must be scalars", call. = FALSE) }

  if (m <= 0 | m == Inf)
  { stop("`m` must be a finite value and greater than 0", call. = FALSE) }

  n = length(y)
  by = mean(y)

  # log-likelihood function
  log_ell = n*( by*log(th) + (m-by)*log(1-th) - log(1-(1-th)^m) )
  el = c(log_ell)
  error = 3

  th_std <- c(th)
  for (k in 1:Maxiter)
  {
    if( error > convergence )
    {
      th = by*(1-(1-th)^m)/m
      th_std <- append(th_std, th)

      log_el = n*( by*log(th) + (m-by)*log(1-th) - log(1-(1-th)^m) )
      el <- append(el, log_el)
      error = abs(el[k+1]-el[k])/(1+abs(el[k]))

      print_err <- error
      print_k <- k
    }
  }
  std_th <- sd(th_std)/sqrt(length(th_std))
  th_t_val <- th/std_th

  ELL = el[length(el)]
  # Rate = ZTB_CRate(m, th)

  #add values of AIC and BIC
  aic <- (2 * length(th)) - (2 * ELL)
  bic <- log(length(y)) * length(th) - 2 * ELL
  info_criteria <- c(AIC=aic, BIC=bic)

  result <- list()
  result$call <- match.call()
  result$print_n <- n
  result$print_k <- print_k
  result$print_err <- print_err
  result$ELL <- ELL
  result$th <- th
  result$std_th <- std_th
  result$th_t_val <- th_t_val
  # result$Rate <- Rate
  result$info_criteria <- info_criteria
  result$convergence <- convergence
  return(result)
}
