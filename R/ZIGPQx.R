#' Title
#'
#' @param phi0 probability
#' @param la scale parameter
#' @param th discrete parameter
#'
#' @return requires a value
#' @export
#'
#' @examples
#' ZIGPQx(0.5,1,0.1)
ZIGPQx <- function(phi0,la,th)
{
  m = length(la)
  df = matrix(0,2*m+1,2*m+1)
  a0 = exp(-sum(la))
  dphi0 = -1/phi0-1/(1-phi0)
  dla = -(1-a0)*(1-phi0)*(1-exp(-la))/la
  dpi = -(1-a0)*(1-phi0)*(1-exp(-la))*la/(th*(1-th))
  dth = c(dphi0,dla,dpi)
  dQ = diag(dth)
  return(dQ)
}
