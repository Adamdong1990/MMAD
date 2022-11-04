summary.ZIGP <- function(object, digits=4, ...){
  cat("Call:\n")
  cat(paste0(deparse(object$call), sep = "\n", collapse = "\n" ))
  cat("\n")
  cat("Model: ")
  cat("Multivariate Compound Zero-Inflated Generalized Poisson Distribution")
  cat("\n\n")
  ZIGP_cnt <- unlist(object$print_n)
  cat("Count of data: ")
  cat(ZIGP_cnt)
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

  m = length(object$la)
  phi0 <- c(round(object$phi0, digits), round(object$std_phi0, digits), round(object$phi0_t_val, digits))
  phi <- c(round(object$phi, digits), round(object$std_phi, digits), round(object$phi_t_val, digits))
  la <- c(round(object$la, digits), round(object$std_la, digits), round(object$la_t_val, digits))
  th <- c(round(object$th, digits), round(object$std_th, digits), round(object$th_t_val, digits))

  phi <- matrix(phi, m, 3)
  la <- matrix(la, m, 3)
  th <- matrix(th, m, 3)

  # rename
  phi_rowname <- c()
  la_rowname <- c()
  th_rowname <- c()
  i = 1
  while (i < m+1)
  {
    p_rowname = paste("phi", i, sep = '_')
    phi_rowname <- append(phi_rowname, p_rowname)
    l_rowname = paste("la", i, sep = '_')
    la_rowname <- append(la_rowname, l_rowname)
    t_rowname = paste("th", i, sep = '_')
    th_rowname <- append(th_rowname, t_rowname)

    i = i + 1
  }

  rownames(phi) <- phi_rowname
  rownames(la) <- la_rowname
  rownames(th) <- th_rowname
  coef_zigp <- rbind(phi0, phi, la, th)
  colnames(coef_zigp) <- c("Estimate", "Std. Error", "t value")
  cat("Coefficients:\n")
  print(coef_zigp)
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
