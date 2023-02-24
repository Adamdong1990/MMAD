summary.GaFrailtyMM <- function(object, digits=4, ...)
{
  cat("Call:\n")
  cat(paste0(deparse(object$call), sep = "\n", collapse = "\n" ))
  cat("\n")
  cat("Model: ")
  cat("Gamma Frailty Survival Models")
  cat("\n\n")

  cat(paste0("n= ", object$N, ", nubetar of events=", object$TotalCen))
  cat("\n")

  error <- object$print_err
  convergence <-  object$convergence
  if (error > convergence)
  {
    cat("Convergence result: Did not converge, try adjust 'Maxiter' or 'convergence'")
    cat("\n")
  }
  else {
    converge_value <- unlist(object$print_err)
    cat("Convergence result: ")
    cat(converge_value)
    cat("\n")
  }

  print_ell <- round(object$ELL, digits)
  cat("\n")

  m = length(object$be)
  theta <- c(round(object$th, digits), round(object$std_th, digits), round(object$ci_th_lower, digits), round(object$ci_th_upper, digits))
  be <- c(round(object$be, digits), round(object$std_be, digits), round(object$ci_be_lower, digits), round(object$ci_be_upper, digits))

  be <- matrix(be, m, 4)

  rownames(be) <- object$namesX

  coef_Frailty <- rbind(theta, be)
  colnames(coef_Frailty) <- c("Estimate", "Std. Error", "Lower 95%-level", "Upper 95%-level")
  cat("Coefficients:\n")
  print(coef_Frailty)
  cat("\n")

  cat("Log Likelihood: ")
  cat(as.character(print_ell))
  cat("\n")
  cat("Optimization Method: ")
  cat("AD technique of MM algorithm")
  cat("\n\n")
}
