#------------------------------------------------------------
# function to estimate stage 1 alpha for exch structure
# for fixed number of categories
#------------------------------------------------------------

#' @export

estalpha1_exch_fixed <- function(mdat, id, Z, nresponse){
  
  match.call()
  
  K <- nresponse
  G1 <- G2 <- 0
  
  for(i in id){
    csizei <- nlevels(as.factor(mdat[mdat$cluster == i,]$unit))
    wi <- 1/csizei
    cveci <- unique(mdat[mdat$cluster == i,]$unit)
    Z_i <- Z[mdat$cluster == i]
    matZ_i <- matrix(Z_i, nrow = K)

    ZZ1 <- 0
    for(j in 1:csizei){
      ZZ1 <- ZZ1 + t(matZ_i[,j]) %*% matZ_i[,j]
    }
    ZZ2 <- sum ( ( colSums(matZ_i) ) ^ 2 )
    
    G1 <- G1 + wi * ZZ1
    G2 <- G2 + wi * ZZ2
  }
  
  thenum <- - ( K - 1) * G1 + sqrt( ( K - 1 ) ^ 2 * G1 ^ 2 - ( G1 - G2 ) * ( ( K - 1 ) * ( G1 * ( K - 1 ) - G2 ) ) )
  theden <- ( K - 1 ) * ( ( K - 1 ) * G1 - G2 )
  
  alpha0 <- thenum / theden

  alpha <- ( alpha0 * ( alpha0 * ( K - 2 ) + 2 ) ) / ( 1 + alpha0 ^ 2 * ( K - 1 ) ) 
  
  R <- toeplitz(c(1, rep(alpha, K - 1)))
  
}

#outb <- estalpha1_exch_fixed(mdat=out[[1]], id =out[[2]], Z=out[[3]], ncategories = out[[4]])
