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
LTNmm <- function(y, a, mu = NULL, sigma = NULL, Maxiter = 2000, convergence = 1e-06, method = "ADMM", ...)
{
  # check 'a'
  if (!is.numeric(a))
  { stop("`a` must be numeric", call. = FALSE) }
  if (length(a) !=1)
  { stop("`a` must be scalars", call. = FALSE) }

  if (a == -Inf | a == Inf)
  { stop("`a` must be a finite value", call. = FALSE) }
  if(any(y < a))
  { stop("observed values below specified truncation `a`", call.=FALSE) }

  # validate and process initial values of mu and sigma
  mu_null <- is.null(mu); sigma_null <- is.null(sigma)
  if (mu_null == TRUE & sigma_null == FALSE)
  { stop("Please give the value of mu", call. = FALSE) }
  if (mu_null == FALSE & sigma_null == TRUE)
  { stop("Please give the value of sigma", call. = FALSE) }

  if(mu_null == TRUE & sigma_null == TRUE)
  {
    lm_mod <- lm(y ~ 1)
    mu <- unname(coef(lm_mod))
    sigma <- log(unname(summary(lm_mod)$sigma))
  }

  error = 3
  n = length(y)
  alpha0 = c(mu,sigma)
  si = sqrt(sigma)
  result <- list()

  # log-likelihood function
  log_ell <- -n*log(sigma)/2-sum((y-mu)^2)/(2*sigma)-n*log(1-pnorm((a-mu)/si))
  el = c(log_ell)

  mu_std <- c(mu)
  sigma_std <- c(sigma)
  for (k in 1:Maxiter)
  {
    if (error > convergence)
    {
      a1 = (a-mu)/si
      w = 1-pnorm(a1)
      s1 = (1-w)/w
      tao = exp(-(a-mu)^2/(2*sigma))/sqrt(2*pi*sigma)
      g = tao/pnorm(a1)
      # deta = sigma - sigma*(a-mu)*g

      mu1 = (mean(y)+s1*(mu-sigma*g))/(1+s1)
      deta = sigma + (mu1-mu)^2 - sigma*(a+mu-2*mu1)*g
      mu = mu1
      sigma = ( sum((y-mu)^2)/n + s1*deta )/(1+s1)
      mu_std <- append(mu_std, mu)
      sigma_std <- append(sigma_std, sigma)

      si = sqrt(sigma)
      log_el <- -n*log(sigma)/2-sum((y-mu)^2)/(2*sigma)-n*log(1-pnorm((a-mu)/si))
      el <- append(el, log_el)
      error = abs(el[k+1]-el[k])/(abs(el[k])+1)

      print_err <- error
      print_k <- k

    }
  }

  std_mu <- sd(mu_std)/sqrt(length(mu_std))
  std_sigma <- sd(sigma_std)/sqrt(length(sigma_std))
  mu_t_val <- mu/std_mu
  sigma_t_val <- sigma/std_sigma

  ELL = el[length(el)]
  alpha = c(mu,sigma)
  Rate = LTN_CRate(a, mu, sigma)

  #add values of AIC and BIC
  aic <- (2 * length(alpha)) - (2 * ELL)
  bic <- log(length(y)) * length(alpha) - 2 * ELL
  info_criteria <- c(AIC=aic, BIC=bic)

  result$call <- match.call()
  result$print_n <- n
  result$print_k <- print_k
  result$print_err <- print_err
  result$ELL <- ELL
  result$mu <- mu
  result$std_mu <- std_mu
  result$mu_t_val <- mu_t_val
  result$sigma <- sigma
  result$std_sigma <- std_sigma
  result$sigma_t_val <- sigma_t_val
  result$Rate <- Rate
  result$info_criteria <- info_criteria
  result$convergence <- convergence
  return(result)
}
