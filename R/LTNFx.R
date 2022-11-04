#' Title
#'
#' @param a left border
#' @param mu mean
#' @param si2 variance
#'
#' @return
#' @export
#'
#' @examples
#' a = 5; mu = 0.5; si2 = 0.5
#' LTNFx(a,mu,si2)
LTNFx = function(a,mu,si2)
{
  a1 = (a-mu)/(sqrt(si2))
  dphi = -a1*exp(-0.5*a1^2)/sqrt(2*pi)
  c = 1-pnorm(a1)
  df = matrix(0,2,2)
  df[1,1] = -1/si2 + dphi/(c*si2) + (dnorm(a1))^2/(si2*c^2)
  df[1,2] = df[2,1] =(a1*dphi -dnorm(a1))/(2*c*sqrt(si2^3)) + a1*(dnorm(a1))^2/(2*c^2*sqrt(si2^3))
  df[2,2] = -1/(2*si2^2) + (a1^2*dphi-a1*dnorm(a1))/(4*c*si2^2)+a1^2*(dnorm(a1))^2/(4*c^2*si2^2)

  return(df)
}
