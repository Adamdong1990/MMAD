#' MM algorithm based on the AD method for Case II interval-censored data
#'
#' @param L Left interval boundary
#' @param R Right interval boundary
#'
#' @return requires a value
#' @export
#'
#' @examples
#' L = c(0.5, 2); R = c(1, 3)
#' IC2mm(L, R)
IC2mm <- function(L, R, Maxiter = 500, convergence = 1e-06, method = "ADMM", ...)
{
  if (length(L) == 0)
  { stop("The length of interval-censored data equals 0", call. = FALSE) }
  if (length(L) != length(R))
  { stop("The amount of interval-censored data (L, R) is not equal", call. = FALSE) }

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

  error = 3
  for (k in 1:Maxiter)
  {
    if (error > convergence)
    {
      A = rowSums(alpha*pp)
      AA = t( matrix(rep(A,each=m),m,n) )
      B = colSums(alpha*pp/AA)
      p = B/sum(B)

      pp = matrix(rep(p,each=n),n,m)
      log_el <- sum( log( rowSums(alpha*pp) ) )
      el <- append(el, log_el)
      error=abs(el[k+1]-el[k])/(1+abs(el[k]))

      print_err <- error
      print_k <- k
    }
  }

  ELL = el[length(el)]

  S = rep(0,m)
  for(i in 1:m)
  {
    S[i] = 1-sum((s<=s[i])*p)
  }

  result <- list()
  result$call <- match.call()
  result$print_n <- n
  result$print_k <- print_k
  result$print_err <- print_err
  result$ELL <- ELL
  result$p <- p
  result$s <- s
  result$S <- S
  result$convergence <- convergence

  return(result)
}
