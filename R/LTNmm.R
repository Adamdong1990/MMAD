#' MM algorithm based on the AD method for Left-truncated normal distribution
#'
#' @param a left border
#' @param mu mean
#' @param si2 variance
#' @param y data
#'
#' @return data
#' @export
#'
#' @examples
#' LTNmm(5, 7, 4, 3)
#' y=c(3, 9); LTNmm(5, 7, 4, y)
LTNmm <- function(a, mu, si2, y)
{
  k = 1
  error = 3
  n = length(y)
  alpha0 = c(mu,si2)
  si = sqrt(si2)
  # log-likelihood function
  log_ell <- -n*log(si2)/2-sum((y-mu)^2)/(2*si2)-n*log(1-pnorm((a-mu)/si))
  el = c(log_ell)

  while (error > 1e-06)
  {
    a1 = (a-mu)/si
    w = 1-pnorm(a1)
    s1 = (1-w)/w
    tao = exp(-(a-mu)^2/(2*si2))/sqrt(2*pi*si2)
    g = tao/pnorm(a1)
    deta = si2 - si2*(a-mu)*g

    mu1 = (mean(y)+s1*(mu-si2*g))/(1+s1)
    deta = si2 + (mu1-mu)^2 - si2*(a+mu-2*mu1)*g
    mu = mu1
    si2 = ( sum((y-mu)^2)/n + s1*deta )/(1+s1)

    si = sqrt(si2)
    log_el <- -n*log(si2)/2-sum((y-mu)^2)/(2*si2)-n*log(1-pnorm((a-mu)/si))
    el <- append(el, log_el)
    error = abs(el[k+1]-el[k])/(abs(el[k])+1)
    k = k+1
  }

  ELL = el[length(el)]
  alpha = c(mu,si2)
  mse = sum( (alpha-alpha0)^2 )/2

  result <- list()
  result$k <- k
  result$ELL <- ELL
  result$alpha <- alpha
  result$mse <- mse
  return(result)
}
