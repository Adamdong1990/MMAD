summary.PET <- function(object, digits=4, ...){
  cat("Call:\n")
  cat(paste0(deparse(object$call), sep = "\n", collapse = "\n" ))
  cat("\n")
  cat("Model: ")
  cat("Poisson Model for Transmission Tomography")
  cat("\n\n")
  PET_cnt <- unlist(object$print_n)
  cat("Count of data: ")
  cat(PET_cnt)
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

  m = length(object$theta)
  theta <- c(round(object$theta, digits), round(object$std_theta, digits), round(object$ci_theta_lower, digits), round(object$ci_theta_upper, digits))

  th <- matrix(theta, m, 4)

  # rename
  theta_rowname <- c()
  i = 1
  while (i < m+1)
  {
    t_rowname = paste("theta", i, sep = '_')
    theta_rowname <- append(theta_rowname, t_rowname)

    i = i + 1
  }

  rownames(th) <- theta_rowname
  coef_pet <- rbind(th)
  colnames(coef_pet) <- c("Estimate", "Std. Error", "Lower 95%-level", "Upper 95%-level")
  cat("Coefficients:\n")
  print(coef_pet)
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
