#' MM algorithm based on the AD method for Zero-truncated binomial distribution(1)
#'
#' @param m number of experiments
#' @param th probability of success
#' @param y
#'
#' @return
#' @export
#'
#' @examples
#' y = c(1, 2)
#' ZTB1mm(5, 0.2, y)
ZTB1mm <- function(m,th,y)
{
  n = length(y)
  by = mean(y)
  th0 = th
  k = 1

  # log-likelihood function
  log_ell = n*( by*log(th) + (m-by)*log(1-th) - log(1-(1-th)^m) )
  el = c(log_ell)
  error = 3

  while( error > 1e-06 )
  {
    th = by*(1-(1-th)^m)/m

    log_el = n*( by*log(th) + (m-by)*log(1-th) - log(1-(1-th)^m) )
    el <- append(el, log_el)
    error = abs(el[k+1]-el[k])/(1+abs(el[k]))
    k = k+1
  }

  ELL = el[length(el)]
  mse = (th-th0)^2

  result <- list()
  result$k <- k
  result$ELL <- ELL
  result$th <- th
  result$mse <- mse
  return(result)
}
