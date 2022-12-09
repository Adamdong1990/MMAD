summary.LTN <- function(object, digits=4, ...){
  cat("Call:\n")
  cat(paste0(deparse(object$call), sep = "\n", collapse = "\n" ))
  cat("\n")
  cat("Model: ")
  cat("Left-Truncated Normal Distribution")
  cat("\n\n")
  LTN_cnt <- unlist(object$print_n)
  cat("Count of data: ")
  cat(LTN_cnt)
  cat("\n")
  converge_cnt <- unlist(object$print_k)
  cat("Number of iterations: ")
  cat(converge_cnt)
  cat("\n")
  error <- object$print_err
  convergence <-  object$convergence
  if (error > convergence)
  {
    cat("Convergence result: Did Not Converge")
    cat("\n")
  }
  else {
    converge_value <- unlist(object$print_err)
    cat("Convergence result: ")
    cat(converge_value)
    cat("\n")
  }

  print_ell <- round(object$ELL, digits)
  print_rate <- round(object$Rate, digits)
  cat("Convergence Rate: ")
  cat(print_rate)
  cat("\n\n")

  mu <- c(round(object$mu, digits), round(object$std_mu, digits), round(object$ci_mu_lower, digits), round(object$ci_mu_upper, digits))
  sigma <- c(round(object$sigma, digits), round(object$std_sigma, digits), round(object$ci_sigma_lower, digits), round(object$ci_sigma_upper, digits))
  coef_ltn <- rbind(mu, sigma)
  colnames(coef_ltn) <- c("Estimate", "Std. Error", "Lower 95%-level", "Upper 95%-level")
  cat("Coefficients:\n")
  print(coef_ltn)
  cat("\n")

  cat("Log Likelihood: ")
  cat(as.character(print_ell))
  cat("\n")
  cat("Information Criterion: ")
  cat(paste0("AIC=", round(object$info_criteria[1], digits), " BIC=", round(object$info_criteria[2], digits)))
  cat("\n")
  cat("Optimization Method: ")
  cat("AD technique of MM algorithm")
  cat("\n\n")
}
