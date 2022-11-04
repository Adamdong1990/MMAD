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
#' ZIGP_CRate(0.5,0.5,1,0.1)
ZIGP_CRate <- function(phi0, phi, la, th)
{
  DF = ZIGPFx(phi0,phi,la,th)
  DGI = solve(ZIGPQx(phi0,phi,la,th))
  A =  DGI%*%DF
  rate = 1-min(eigen(A)$values)
  return(rate)
}
