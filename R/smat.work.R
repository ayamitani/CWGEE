# ------------------------------------------------------------
# function to construct S matrix
# ------------------------------------------------------------

smat.work <- function(coef,categs){
  
  # matrix S
  categs1 <- categs-1
  S_mat <- matrix(nrow=categs1,ncol=categs1)
  for (j in 1:categs1){
    for (k in 1:categs1){
      if (j>=k){
        S_mat[j,k] <- sqrt(exp(coef[k]-coef[j]))
      }
      if (j<k){
        S_mat[j,k] <- sqrt(exp(coef[j]-coef[k]))
      }
    } # k
  } # j
  
  # output
  list(S_mat=S_mat)
}