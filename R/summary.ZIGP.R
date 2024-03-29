#' Summary of parameter estimates of a Type I multivariate ZIGP model
#'
#' @description This function returns the result of the \code{ZIGPMM} function
#'
#'
#' @aliases summary.ZIGP
#' @usage \method{summary}{ZIGP}(object, digits = 4, ...)
#' @param object Output from a call to ZIGP.
#' @param digits The desired number of digits after the decimal point. Default of 4 digits is used.
#' @param ... Additional arguments
#'
#' @return Summary for \code{ZIGPMM} objects.
#' @seealso \code{\link{ZIGPMM}}
#' @keywords methods
#' @method summary ZIGP
#' @export
#'
#' @examples
#'
#'
#' x1 <- c(0, 0, 0,38, 0,19,25, 0,25, 0)
#' x2 <- c(0, 0, 0,23, 0,51,24, 0,10, 0)
#' y <- cbind(x1, x2)
#' phi0 = 0.5; la = rep(1,2); th = rep(0.1,2)
#' result <- ZIGPMM(y, phi0, la, th)
#'
#' summary(result,digits=4)
#'
#'
#'
summary.ZIGP <- function(object, digits = 4, ...) {
  cat("Call:\n")
  cat(paste0(deparse(object$call), sep = "\n", collapse = "\n"))
  cat("\n")
  cat("Model: ")
  cat("Type I Multivariate Zero-Inflated Generalized Poisson Distribution")
  cat("\n\n")
  ZIGP_cnt <- unlist(object$print_n)
  cat("Count of data: ")
  cat(ZIGP_cnt)
  cat("\n")
  converge_cnt <- unlist(object$print_k)
  cat("Number of iterations: ")
  cat(converge_cnt)
  cat("\n")

  print_ell <- round(object$ELL, digits)
  print_rate <- round(object$Rate, digits)
  cat("Convergence Rate: ")
  cat(print_rate)
  cat("\n\n")

  m <- length(object$la)
  phi0 <- c(round(object$phi0, digits), round(object$std_phi0, digits),
            round(object$ci_phi0_lower, digits), round(object$ci_phi0_upper,
                                                       digits))
  la <- c(round(object$la, digits), round(object$std_la, digits), round(object$ci_la_lower,
                                                                        digits), round(object$ci_la_upper, digits))
  th <- c(round(object$th, digits), round(object$std_th, digits), round(object$ci_th_lower,
                                                                        digits), round(object$ci_th_upper, digits))

  la <- matrix(la, m, 4)
  th <- matrix(th, m, 4)

  # rename
  la_rowname <- c()
  th_rowname <- c()
  i <- 1
  while (i < m + 1) {
    l_rowname <- paste("la", i, sep = "_")
    la_rowname <- append(la_rowname, l_rowname)
    t_rowname <- paste("th", i, sep = "_")
    th_rowname <- append(th_rowname, t_rowname)

    i <- i + 1
  }

  rownames(la) <- la_rowname
  rownames(th) <- th_rowname
  coef_zigp <- rbind(phi0, la, th)
  colnames(coef_zigp) <- c("Estimate", "Std. Error", "Lower 95%-level",
                           "Upper 95%-level")
  cat("Coefficients:\n")
  print(coef_zigp)
  cat("\n")

  cat("Log Likelihood: ")
  cat(as.character(print_ell))
  cat("\n")
  cat("Information Criterion: ")
  cat(paste0("AIC=", round(object$info_criteria[1], digits), " BIC=",
             round(object$info_criteria[2], digits)))
  cat("\n")
  cat("Optimization Method: ")
  cat("AD technique of MM algorithm")
  cat("\n\n")
}
