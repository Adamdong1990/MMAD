#' Title
#'
#' @param a left border
#' @param mu mean
#' @param si2 variance
#'
#' @return requires a value
#' @export
#'
#' @examples
#' a = 5; mu = 0.5; si2 = 0.5
#' LTNQx(a, mu, si2)
LTNQx = function(a, mu, si2)
{
  a1 = (a-mu)/(sqrt(si2))
  c = 1-pnorm(a1)
  s1 = (1-c)/c
  tao = exp(-(a-mu)^2/(2*si2))/sqrt(2*pi*si2)
  g = tao/pnorm(a1)

  dQ = matrix(0,2,2)
  dQ[1,1] = -(1+s1)/si2
  dQ[1,2] = dQ[2,1] = s1*g/si2 - dnorm(a1)/(c*sqrt(si2^3))
  dQ[2,2] = -(1+s1)/(2*si2^2) + s1*(a-mu)*g/(si2^2) - a1*dnorm(a1)/(c*si2^2)

  return(dQ)
}
