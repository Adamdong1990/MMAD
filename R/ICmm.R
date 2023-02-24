ICmm <- function(L, R, Lin=NULL, Rin=NULL, Maxiter = 1000, convergence = 1e-06, method = "ADMM", ...)
{
  if (length(L) == 0)
  { stop("The length of interval-censored data equals 0", call. = FALSE) }
  if (length(L) != length(R))
  { stop("The amount of interval-censored data (L, R) is not equal", call. = FALSE) }

  AI<-Aintmap(L,R,Lin,Rin)
  alpha<-AI$A
  intmap<-AI$intmap

  if (any(apply(alpha,1,sum)==0)) stop(" `alpha` row all zeros. Appears that there are some R<L")
  n<-dim(alpha)[[1]]
  m<-dim(alpha)[[2]]
  alpha<-matrix(alpha,n,m)
  if (m==1)
  {
    p <- 1
    print_err <- 0
    print_k <- 1
    ELL <- sum( log( rowSums(alpha*pp) ) )
    ### fix error 1/24/2011: needed to change name from numit to count in emout list
    # emout<-list(error=0,count=0,converge=TRUE,message="normal convergence")
    # anypzero<-FALSE
  }
  else
  {
    p = rep(1/m,m)
    pp =  matrix(rep(p,each=n),n,m)
    log_ell = sum( log( rowSums(alpha*pp) ) )
    el = c(log_ell)

    error = 3
    for (k in 1:Maxiter)
    {
      if (error > convergence)
      {
        A1 = rowSums(alpha*pp)
        AA = t( matrix(rep(A1,each=m),m,n) )
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
  }

  ps = c(0,p)
  s = c(0,AI$intmapR)
  S = rep(0,m+1)
  for(i in 1:(m+1))
  {
    S[i] = 1-sum((s<=s[i])*ps)
  }

  result <- list()
  result$call <- match.call()
  result$intmap <- intmap
  result$print_k <- print_k
  result$print_err <- print_err
  result$ELL <- ELL
  result$p <- p
  result$s <- s
  result$S <- S
  result$convergence <- convergence

  return(result)
}
