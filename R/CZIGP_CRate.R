#' Title
#'
#' @param phi0 probability
#' @param phi probability
#' @param la scale parameter
#' @param th discrete parameter
#'
#' @return requires a value
#' @export
#'
#' @examples
#' CZIGP_CRate(0.5,0.5,1,0.1)
CZIGP_CRate <- function(phi0, phi, la, th)
{
  DF = CZIGPFx(phi0,phi,la,th)
  DGI = solve(CZIGPQx(phi0,phi,la,th))
  A =  DGI%*%DF
  rate = 1-min(eigen(A)$values)
  return(rate)
}
