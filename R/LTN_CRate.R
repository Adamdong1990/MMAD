#' Title
#'
#' @param a left border
#' @param mu mean
#' @param sig2 variance
#'
#' @return requires a value
#' @export
#'
#' @examples
#' a = 5; mu = 0.5; si2 = 0.5
#' LTN_CRate(a, mu, si2)
LTN_CRate <- function(a, mu, sig2)
{
  DF = LTNFx(a,mu,si2)
  DGI = solve(LTNQx(a,mu,si2))
  A =  DGI%*%DF
  rate = 1-min(eigen(A)$values)
  return(rate)
}
