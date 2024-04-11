#' Generate sample from Multivariate Zero-inflated Compound Generalized Poisson Distribution

#' @param phi0 Probability of taking value 0 for outer Bernoulli distribution
#' @param phi Vector of probability of taking value 0 for inner Bernoulli distributions
#' @param la Lambda for Poison distribution
#' @param th Theta for Poison distribution
#' @param n Number of samples generated
#'
#' @return Vector of samples with dimension n
#'
#' @export

CZIGP_sample <- function(phi0, phi, la, th, n) {
  GP <- function(x, la0, th0) {
    b <- gamma(x + 1)
    if (th0 < 1 && th0 >= 0) {
      p <- exp(-la0 - th0 * x) * la0 * (la0 + th0 * x)^(x - 1)/b
    }
    if (th0 < 0 && th0 > -1 && la0 + th0 * x <= 0) {
      p <- 0
    }
    if (th0 < 0 && th0 > -1 && la0 + th0 * x > 0) {
      p <- exp(-la0 - th0 * x) * la0 * (la0 + th0 * x)^(x - 1)/b
    }
    return(p)
  }

  m <- length(la)
  Y <- matrix(0, n, m)
  P <- matrix(0, 100, m)
  for (j in 1:m) {
    for (i in 0:99) {
      P[i + 1, j] <- GP(i, la[j], th[j])
    }
    a <- 0:99
    X <- sample(a, n, P[, j], replace = TRUE)
    Z <- sample(c(0, 1), n, c(phi[j], 1 - phi[j]), replace = TRUE)
    Y[, j] <- X * Z
  }
  z0 <- sample(c(0, 1), n, c(phi0, 1 - phi0), replace = TRUE)
  zz <- matrix(z0, n, m)
  y <- zz * Y
  return(y)
}

