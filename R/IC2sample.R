#' Generate sample for Type II interval-censored data

#' @param la Hazard rate for observed data
#' @param n Number of samples generated
#' @return List of n objects containing left and right bounds for observed data
#'
#' @export

IC2sample <- function(la, n) {

  u <- runif(n, 0, 1)
  w <- runif(n, 0, 1)
  U <- -log(u)/la
  W <- -log(w)/la
  L <- pmin(U, W)
  R <- pmax(U, W)

  return(list(L = L, R = R))
}
