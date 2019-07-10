#------------------------------------------------------------
# function to estimate stage 2 alpha for ar1 structure
#------------------------------------------------------------

estalpha2_ar1 <- function(alpha0){
  alpha <- as.numeric( 2 * alpha0 / ( 1 + alpha0 ^ 2 ) )
  return(alpha)
}