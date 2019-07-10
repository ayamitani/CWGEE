#----------------------------------
# Print output
#----------------------------------

#' @export

print.cwgee <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}