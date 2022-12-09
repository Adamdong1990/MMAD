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
#' ZIGP_Rate(0.5,1,0.1)
ZIGP_Rate <- function(phi0, la, th)
{
  DF = ZIGPFx(phi0,la,th)
  DGI = solve(ZIGPQx(phi0,la,th))
  A =  DGI%*%DF
  rate = 1-min(eigen(A)$values)
  return(rate)
}
