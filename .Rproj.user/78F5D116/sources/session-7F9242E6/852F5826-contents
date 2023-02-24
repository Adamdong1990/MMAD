summary.cox <- function(object, digits=4, ...){
  cat("Call:\n")
  cat(paste0(deparse(object$call), sep = "\n", collapse = "\n" ))
  cat("\n")
  cat("Model: ")
  cat("Cox Model with Right-censored")
  cat("\n\n")

  cat(paste0("n= ", object$n, ", number of events=", object$event))
  cat("\n\n")

  m = length(object$be)
  be = round(object$be, digits)
  se = round(object$se, digits)
  ci_normal = c( round(object$ci_lower, digits), round(object$ci_upper, digits) )
  beta <- cbind(matrix(be, ncol=1), matrix(se, ncol=1), matrix(ci_normal, ncol=2))

  # non-normal
  beta_non <- cbind( matrix(be, ncol=1), matrix(se, ncol=1), round(object$ci, digits) )

  # rename
  beta_rowname <- c()
  i = 1
  while (i < m+1)
  {
    b_rowname = paste("beta", i, sep = '_')
    beta_rowname <- append(beta_rowname, b_rowname)

    i = i + 1
  }

  rownames(beta) <- beta_rowname
  rownames(beta_non) <- beta_rowname
  coef <- rbind(beta)
  coef_non <- rbind(beta_non)
  colnames(coef) <- c("Estimate", "Std. Error", "Lower 95%-level", "Upper 95%-level")
  cat("Coefficients:\n")
  cat("bootstrap replications are normally distributed:\n")
  print(coef)
  cat("\n")
  cat("bootstrap replications are non-normally distributed:\n")
  colnames(coef_non) <- c("Estimate", "Std. Error", "Lower 95%-level", "Upper 95%-level")
  print(coef_non)
  cat("\n")

  method <- ifelse(object$method=="Pro", "Profile", "Non Profile" )
  cat("Optimization Method: ")
  cat("AD technique of MM algorithm(", method, ")")
  cat("\n\n")
}
