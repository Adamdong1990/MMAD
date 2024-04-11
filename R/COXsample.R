#' Generate sample from COX model

#' @param n Number of samples generated
#' @param be Coefficients for COX regression
#' @param la Baseline lambda for COX regression (we take constant value here)
#' @param n Number of samples generated
#'
#' @return List with length n including covariate x, censoring indicator I and event time t
#'
#' @export


COXsample <- function(n, be, la) {
  q <- length(be)
  x <- matrix(rnorm(n * q, -1, 1), n, q)
  u <- runif(n)
  t <- -log(u)/(la * exp(x %*% be))
  cen <- 3.7
  I <- 1 * (t <= cen)
  t <- pmin(t, cen)
  return(list(x = x, I = I, t = t))
}
