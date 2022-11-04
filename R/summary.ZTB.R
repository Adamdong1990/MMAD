summary.ZTB <- function(object, digits=4, ...){
  cat("Call:\n")
  cat(paste0(deparse(object$call), sep = "\n", collapse = "\n" ))
  cat("\n")
  cat("Model: ")
  cat("Zero-Truncated Binomial Distribution")
  cat("\n\n")
  ZTB_cnt <- unlist(object$print_n)
  cat("Count of data: ")
  cat(ZTB_cnt)
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
  # print_rate <- round(object$Rate, digits)
  # cat("Convergence Rate: ")
  # cat(print_rate)
  cat("\n\n")

  th <- c(round(object$th, digits), round(object$std_th, digits), round(object$th_t_val, digits))
  coef_ltn <- matrix(th, 1, 3)
  rownames(coef_ltn) <- "th "
  colnames(coef_ltn) <- c("Estimate", "Std. Error", "t value")
  cat("Coefficients:\n")
  print(coef_ltn)
  cat("\n")

  # cat("MSE: ")
  # cat(round(object$mse, digits))
  # cat("\n")
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
