#' Title
#'
#' @param data a data frame containing the variables of the model
#' @param time for right censored data
#' @param status The status indicator, normally 0=alive, 1=dead.
#' @param be a vector of regression parameters
#' @param Maxiter The maximum iterations
#' @param convergence convergence value
#' @param bootstrap bootstrap
#' @param b number of bootstrap iterations for the estimation of the standard errors
#' @param method Profile or Non Profile
#' @param ...
#'
#' @return requires a value
#' @export
#'
#' @examples
#' data = matrix(rnorm(10*1,-1,1),10,1)
#' be = 0.5
#' u=runif(10); time=-log(u)/(la*exp(x%*%be))
#' cen=3.7; status=1*(t<=cen)
#' coxmm(data, time, status, be)
coxmm <- function(data, time, status, be, Maxiter = 500, convergence = 1e-06, bootstrap = TRUE, b = 100, method = "Pro", ...)
{
  o.data <- data
  o.time <- time
  o.status <- status
  o.be <- be
  se <- ci <- NULL

  n=dim(data)[1]
  event=sum(status)
  la=rep(1/n,n)
  eps=1

  for (k in 1:Maxiter)
  {
    if (eps > convergence)
    {
      if(method == "Pro")
      { re=Profile(data, status, time, be, la) }
      else if (method == "NPro")
      { re=Non_Profile(data, status, time, be, la) }

      be=re$newbe
      la=re$newla
      LA=re$newLA
      eps=re$eps
    }
  }

  if (bootstrap == TRUE)
  {
    result_se <- std.est(
      data = o.data,
      time = o.time,
      status = o.status,
      be = o.be,
      Maxiter = Maxiter,
      convergence = convergence,
      bootstrap = bootstrap,
      b = b,
      method = method
    )

    se <- result_se$se
    ci <- result_se$ci
  }

  # construct a Wald-type bootstrap confidence interval, if approximately normally distributed
  # Wald-type bootstrap CI is
  ci_lower <- be - se*1.96
  ci_upper <- be + se*1.96

  result <- list()
  result$call <- match.call()
  result$n <- n
  result$event <- event
  result$be <- be
  result$la <- la
  result$LA <- LA
  result$eps <- eps
  result$convergence <- convergence
  result$se <- se
  result$ci <- ci
  result$ci_lower <- ci_lower
  result$ci_upper <- ci_upper
  result$method <- method

  return(result)
}
