summary.IC2 <- function(object, digits=4, ...)
{
  cat("Call:\n")
  cat(paste0(deparse(object$call), sep = "\n", collapse = "\n" ))
  cat("\n")
  cat("Model: ")
  cat("Case II Interval-Censored Data")
  cat("\n\n")
  IC2_cnt <- unlist(object$print_n)
  cat("Count of data: ")
  cat(IC2_cnt)
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

  m <- length(object$p)
  p <- c(round(object$p, digits), round(object$std_p, digits), round(object$p_t_val, digits))
  p <- matrix(p, m, 3)

  # rename
  p_rowname <- c()
  i = 1
  while (i < m+1)
  {
    p_name = paste("p", i, sep = '_')
    p_rowname <- append(p_rowname, p_name)
    i = i + 1
  }
  rownames(p) <- p_rowname
  colnames(p) <- c("Estimate", "Std. Error", "t value")
  cat("Coefficients:\n")
  print(p)
  cat("\n")

  cnt_tbl <- c(round(object$MAE, digits), round(object$MSE, digits))
  names(cnt_tbl) <- c(" MAE", " MSE")
  print(cnt_tbl)
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
