#' MM algorithm based on the AD method for Case II interval-censored data
#'
#' @param L Left interval boundary
#' @param R Right interval boundary
#'
#' @return
#' @export
#'
#' @examples
#' L = c(0.5, 2); R = c(1, 3)
#' IC2mm(L, R)
IC2mm <- function(L, R)
{
  if (length(L) == 0)
    stop("The length of interval-censored data equals 0")
  if (length(L) != length(R))
    stop("The amount of interval-censored data (L, R) is not equal")

  LR = c(L,R)
  ss = sort(LR)
  s = unique(ss)
  m = length(s)
  n = length(L)
  alpha = matrix(0,n,m)
  for(i in 1:n)
  {
    alpha[i,] = (s>L[i])*(s<=R[i])
  }

  p = rep(1/m,m)
  pp =  matrix(rep(p,each=n),n,m)
  log_ell = sum( log( rowSums(alpha*pp) ) )
  el = c(log_ell)

  k=1
  error = 3
  while( error > 1e-06 )
  {
    A = rowSums(alpha*pp)
    AA = t( matrix(rep(A,each=m),m,n) )
    B = colSums(alpha*pp/AA)
    p = B/sum(B)

    pp = matrix(rep(p,each=n),n,m)
    log_el <- sum( log( rowSums(alpha*pp) ) )
    el <- append(el, log_el)
    error=abs(el[k+1]-el[k])/(1+abs(el[k]))
    k = k+1
  }

  ELL = el[length(el)]
  S = rep(0,m)
  for(i in 1:m)
  {
    S[i] = 1-sum((s<=s[i])*p)
  }

  M =  max(abs(S-exp(-0.5*s)))
  t1 = c(0,s[1:(m-1)])
  t = s-t1
  AM = sum( (S-exp(-0.5*s))^2*t )

  result <- list()
  result$k <- k
  result$ELL <- ELL
  result$p <- p
  result$M <- M
  result$AM <- AM

  return(result)
}
