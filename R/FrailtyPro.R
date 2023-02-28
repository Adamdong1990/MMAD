#' Fitting Gamma frailty models with the profile MM algorithm
#'
#' @param N Total amount of observations.
#' @param q The number of unknown regression parameters.
#' @param x A vector of covariates.
#' @param d Censored indicator: 0 or 1, default 0= censored.
#' @param a The number of clusters.
#' @param b The number of members in each cluster.
#' @param lambda Baseline hazard rate.
#' @param theta The Variance of frailty factors subject to Gamma distribution.
#' @param beta A vector of unknown regression parameters.
#' @param vy The survival time in vector form, e.g. \code{as.vector(time)}.
#' @param vd The Censored indicator in vector form, e.g. \code{as.vector(status)}.
#'
#' @return The final estimate of theta, beta and lambda.
#' @export
#'
#' @examples
#' a=6; b=2; q=2; N=a*b
#' x <- array(NA, c(a,b,q))
#' x1 <- c(28,48,32,31,10,16, 28,48,32,32,10,17)
#' x2 <- c(1,2,1,2,1,2 ,1,2,1,2,1,2)
#' x[,,1] <- matrix(x1, 6, 2)
#' x[,,2] <- matrix(x2, 6, 2)
#' y1 <- c(8, 23, 22,447,30,24,16,13,28,318,12,245)
#' d1 <- c(1,2,1,2,1,2 ,1,2,1,2,1,2)
#' y <- matrix(y1, 6, 2)
#' d <- matrix(d2, 6, 2)
#' vy <- as.vector(y)
#' vd <- as.vector(d)
#' beta = rep(1, q); theta = 1; lambda = rep(1/N, N)
#' FrailtyPro(N, q, x, d, a, b, lambda, theta, beta, vy, vd)
#'
FrailtyPro <- function(N, q, x, d, a, b, lambda, theta, beta, vy, vd)
{
  La = ( cumsum( lambda[order(vy)] ) )[rank(vy)]
  La = matrix(La,a,b)
  A = 1/theta+rowSums(d)
  BE = array(rep(beta,each=a*b),c(a,b,length(beta)))
  C = 1/theta+rowSums(La*exp(apply(x*BE,c(1,2),sum)))
  AC = matrix(A/C,a,b)
  E_0 = as.vector( AC*exp(apply(x*BE,c(1,2),sum)) )
  SUM_0 = cumsum( (E_0[order(vy)])[seq(N,1,-1)] )
  SUM_0 = ( SUM_0[seq(N,1,-1)] )[rank(vy)]

  # beta #
  for(p in 1:q)
  {
    E_1 = as.vector( AC*x[,,p]*exp(apply(x*BE,c(1,2),sum)) )
    AVE_X = apply(abs(x),c(1,2),sum)/abs(x[,,p])
    E_2 = as.vector( AC*AVE_X*x[,,p]^2*exp(apply(x*BE,c(1,2),sum))  )

    SUM_1 = cumsum( (E_1[order(vy)])[seq(N,1,-1)] )
    SUM_1 = ( SUM_1[seq(N,1,-1)] )[rank(vy)]
    SUM_2 = cumsum( (E_2[order(vy)])[seq(N,1,-1)] )
    SUM_2 = ( SUM_2[seq(N,1,-1)] )[rank(vy)]

    DE_1 = sum(d*x[,,p]) - sum(vd*SUM_1/SUM_0)
    DE_2 = -sum(vd*SUM_2/SUM_0)

    beta[p] = beta[p] - DE_1/DE_2
  }

  # th #
  Q01 = a*(digamma(1/theta)+log(theta)-1)/(theta^2) + sum(A/C-digamma(A)+log(C))/(theta^2)
  Q02 = a*(3-2*digamma(1/theta)-2*log(theta))/(theta^3)+2*sum(digamma(A)-log(C)-A/C)/(theta^3) - a*trigamma(1/theta)/(theta^4)
  th = theta - Q01/Q02
  if(th>0) { theta = th }

  # la #
  A = 1/theta +rowSums(d)
  BE = array(rep(beta,each=a*b),c(a,b,length(beta)))
  C = 1/theta+rowSums(La*exp(apply(x*BE,c(1,2),sum)))
  AC = matrix(A/C,a,b)
  E_0 = as.vector( AC*exp(apply(x*BE,c(1,2),sum)) )
  SUM_0 = cumsum( (E_0[order(vy)])[seq(N,1,-1)] )
  SUM_0 = ( SUM_0[seq(N,1,-1)] )[rank(vy)]
  lambda = vd/SUM_0
  TB = c(theta,beta)
  return(list(TB=TB,lambda=lambda))
}
