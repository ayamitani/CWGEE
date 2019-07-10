#------------------------------------------------------------
# function to construct C matrix
#------------------------------------------------------------

corr.str <- function(d, time.str, rho){
  
  match.call()
  types <- c("ind", "exch", "ar1")
  if (which(types == time.str) == 1){
    corr.mat <- matrix(0, ncol = d, nrow = d)
    diag(corr.mat) <- 1
  }
  if (which(types == time.str) == 2)
    corr.mat <- toeplitz(c(1, rep(rho, d - 1)))
  if (which(types == time.str) == 3)
    corr.mat <- suppressWarnings(matrix( rep( rho^abs(-d:d), d )[-(1:d)], ncol=d )[1:d,1:d])
  
  return(corr.mat)
  
}