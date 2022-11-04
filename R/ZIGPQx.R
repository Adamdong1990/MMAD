#' Title
#'
#' @param phi0 probability
#' @param phi probability
#' @param la scale parameter
#' @param th discrete parameter
#'
#' @return
#' @export
#'
#' @examples
#' ZIGPQx(0.5,0.5,1,0.1)
ZIGPQx = function(phi0, phi, la, th)
{
  m = length(phi)
  df = matrix(0,3*m+1,3*m+1)
  r1 = phi0+(1-phi0)*prod(phi+(1-phi)*exp(-la))
  be = phi+(1-phi)*exp(-la)
  dphi0 = -1/phi0-1/(1-phi0)
  dphi = -(r1-phi0)/(be*phi)-(r1-phi0)*exp(-la)/(be*(1-phi))-(1-r1)/phi-(1-r1)/(1-phi)
  dla = -(1-r1)*(1-phi)*(1-exp(-la))/la
  dth = -(1-r1)*(1-phi)*(1-exp(-la))*la/(th*(1-th))
  dqth = c(dphi0,dphi,dla,dth)
  dQ = diag(dqth)
  return(dQ)
}
