#' The MM algorithm based on AD technology for gamma frailty model
#'
#' @param formula A formula object, which contains on the left hand side an object of the type \code{Surv}
#' and on the right hand side a \code{+cluster(id)} statement.
#' e.g. \code{formula=Surv(time, status) ~ x + cluster(id)}
#' @param data A \code{data.frame} in which to interpret the variables named in the formula.
#' @param beta A vector of unknown regression parameters, default is \code{NULL}.
#' If is \code{NULL}, then make all \code{beta=1} during calculation.
#' @param theta The Variance of frailty factors subject to Gamma distribution, default is \code{NULL}.
#' If is \code{NULL}, then let \code{theta=1} during calculation.
#' @param lambda Baseline hazard rate, default set to \code{NULL}. If is \code{NULL},
#' then let each \code{lambda} equals to \code{1/N} during calculation, which \code{N} is the number of observed.
#' @param Maxiter The maximum number of iterations; default=1000.
#' @param convergence Specify the convergence criterion, the default is 1e-6.
#' @param ... Additional arguments
#'
#' @return An object of class \code{GaFrailtyMM} that contains the following fields:
#' \item{N}{list of total observations}
#' \item{TotalCen}{list of censored observations}
#' \item{print_k}{The iterations}
#' \item{print_err}{Convergence result}
#' \item{ELL}{The Log Likelihood value}
#' \item{th}{\eqn{\theta}}
#' \item{std_th}{The standard deviation of the estimated \eqn{\theta}}
#' \item{ci_th_lower, ci_th_upper}{The likelihood-based 95\% confidence interval for the \eqn{\theta}}
#' \item{be}{\eqn{\beta}}
#' \item{std_be}{The standard deviation of the estimated \eqn{\beta}}
#' \item{ci_be_lower, ci_be_upper}{The likelihood-based 95\% confidence interval for the \eqn{\beta}}
#' \item{convergence}{The final convergence value}
#' \item{namesX}{The variable name}
#' @export
#'
#' @examples
#' gf = GaFrailtyMM(Surv(time, status) ~ age + sex + cluster(id), kidney)
#' summary.GaFrailtyMM(gf)
GaFrailtyMM <- function(formula, data, beta=NULL, theta=NULL, lambda=NULL, Maxiter = 1000, convergence = 1e-06, ...)
{
  if (missing(formula) || class(formula) != "formula")
  { stop("'formula' is missing or incorrect") }
  if (missing(data) || class(data) != "data.frame")
  { stop("'data' is missing or not an object of type data.frame") }

  result <- list()
  mf <- model.frame(formula, data)
  m <- model.matrix(formula, data)

  # Identify the cluster
  cluster_id <- grep("cluster", names(mf))
  if(length(cluster_id) != 1)
  { stop("misspecified or non-specified cluster") }
  id <- mf[[cluster_id]]
  mxid <- table(id)
  uid <- unique(as.numeric(mxid))
  if (length(uid) != 1)
  { stop("Unequal nubetar of events per cluster") }
  a = length( names(mxid) )
  b = uid
  N = a*b

  cluster_X <- grep("cluster", colnames(m))
  X <- m[,-c(1, cluster_X), drop=FALSE]
  namesX <- colnames(X)
  q = length(namesX)

  nord = order(id)
  idx = seq(1,N, b)
  x <- array(NA, c(a,b,q))

  for (i in 1:q)
  {
    for (j in 1:b)
    {
      idx = seq(j,N, b)
      x[,,i][,j] = c(X[,i][idx])
    }
  }

  y = matrix(mf[[1]][nord, 1], c(a, b), byrow = TRUE)
  d = matrix(mf[[1]][nord, 2], c(a, b), byrow = TRUE)

  TotalCen = sum(d)
  vy = as.vector(y)
  vd = as.vector(d)

  if (is.null(beta) == TRUE)
  { beta = rep(1, q) }
  if (is.null(theta) == TRUE)
  { theta = 1 }
  if (is.null(lambda) == TRUE)
  { lambda = rep(1/N, N) }

  log_ell = LogLik(lambda, x, d, vy, vd, theta, beta, a, b)
  ell = c(log_ell)

  error = 3

  for (k in 1:Maxiter)
  {
    if (error > convergence)
    {
      Fa = FrailtyPro(N, q, x, d, a, b, lambda, theta, beta, vy, vd)
      lambda = Fa$lambda
      TB = Fa$TB
      r_n = TB[2:(q+1)]-beta
      if(TB[1]>0) theta=TB[1]

      FFa = FrailtyPro(N, q, x, d, a, b, lambda, theta, TB[2:(q+1)], vy, vd)
      lambda = FFa$lambda
      TB1 = FFa$TB
      v_n = TB1[2:(q+1)]-2*TB[2:(q+1)]+beta
      al = (sum(r_n*v_n))/(sum(v_n^2))
      if(al>-1) al=-1

      beta=beta-2*al*r_n+al^2*v_n

      if(TB1[1]>0) theta = TB1[1]

      La = ( cumsum( lambda[order(vy)] ) )[rank(vy)]
      La = matrix(La,a,b)
      A = 1/theta+rowSums(d)
      BE = array(rep(beta,each=a*b),c(a,b,length(beta)))
      C = 1/theta+rowSums(La*exp(apply(x*BE,c(1,2),sum)))
      AC = matrix(A/C,a,b)
      E_0 = as.vector( AC*exp(apply(x*BE,c(1,2),sum)) )
      SUM_0 = cumsum( (E_0[order(vy)])[seq(N,1,-1)] )
      SUM_0 = ( SUM_0[seq(N,1,-1)] )[rank(vy)]
      lambda = vd/SUM_0

      log_el=LogLik(lambda, x, d, vy, vd, theta, beta, a, b)
      ell <- append(ell, log_el)
      error=abs(ell[k+1]-ell[k])/(1+abs(ell[k]))

      print_err <- error
      print_k <- k
    }

  }

  ELL = ell[length(ell)]
  Hessian_Matrix <- HessianF(N, q, x, d, vd, a, b, lambda, theta, beta, vy)
  std_err <- sqrt(diag( solve(Hessian_Matrix) ))
  std_th <- std_err[1]
  std_be <- std_err[2: (1+q)]

  # confidence intervals
  ci_th_lower <- theta - 1.96*std_th
  ci_th_upper <- theta + 1.96*std_th

  ci_be_lower <- beta - 1.96*std_be
  ci_be_upper <- beta + 1.96*std_be

  result$call <- match.call()
  result$N <- N
  result$TotalCen <- TotalCen
  result$print_k <- print_k
  result$print_err <- print_err
  result$ELL <- ELL
  result$th <- theta
  result$std_th <- std_th
  result$ci_th_lower <- ci_th_lower
  result$ci_th_upper <- ci_th_upper
  result$be <- beta
  result$std_be <- std_be
  result$ci_be_lower <- ci_be_lower
  result$ci_be_upper <- ci_be_upper
  result$convergence <- convergence
  result$namesX <- namesX
  return(result)
}
