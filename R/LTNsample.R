#' Generate sample from Left Truncated Normal Distribution

#' @param a Left truncation point
#' @param mu Distribution mean for normal distribution before left truncation
#' @param si2 Distribution variance for normal distribution before left truncation
#' @param n Number of samples generated
#'
#' @return Vector of samples with dimension n
#'
#' @export

LTN_sample <- function(a, mu, si2, n) {
  x <- rep(0, n)
  for (i in 1:n) {
    repeat {
      z <- rnorm(1, mu, sqrt(si2))
      if (z >= a) {
        break
      }
    }
    x[i] <- z
  }
  return(x)
}
