#----------------------------------
# Print summary
#----------------------------------

#' @export

print.summary.cwgee <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients)
  cat("\nWorking correlation structure:\n")
  print(x$corr.str)
  cat("\nCorrelation matrix estimate:\n")
  print(x$corr.matrix)
  cat("\nTemporal working correlation structure:\n")
  print(x$time.str)
  cat("\nCorrelation estimate:\n")
  print(x$corr.coef)
  cat("\nNumber of iterations:\n")
  print(x$niter)
}