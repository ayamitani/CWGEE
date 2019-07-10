#------------------------------------------------------------
# function to estimate stage 1 alpha for ar1 structure
#------------------------------------------------------------


estalpha1_ar1 <- function(mdat, id, Z, Sinv, ncategories1){
  
  Fa <- Fb <- Fc <- 0
  for (i in id){
    csizei <- nlevels(as.factor(mdat[mdat$id == i,]$cluster.var))
    wi <- 1/csizei
    cveci <- unique(mdat[mdat$id == i,]$cluster.var)
    S1_j <- S2_j <- S2_ja <- S2_jb <- 0
    for (j in cveci){
      t_ij <- nlevels(as.factor(mdat[mdat$id == i & mdat$cluster.var == j,]$time.var))
      wwij <- 1 / t_ij
      Z_ij <- Z[mdat$id == i & mdat$cluster.var == j]
      matZ_ij <- matrix(Z_ij, nrow = ncategories1)
      if (t_ij > 1){
        for(c in 1:(t_ij-1)){
          S1_j <- S1_j + wwij * ( t(matZ_ij[,c]) %*% Sinv %*% matZ_ij[,c+1] )
        }
        if (t_ij == 2){
          for(c in 1:t_ij){
            S2_j <- S2_j + wwij * ( t(matZ_ij[,c]) %*% Sinv %*% matZ_ij[,c] )
          }
        }else{
          for(c in 1:t_ij){
            S2_ja <- S2_ja + wwij * ( t(matZ_ij[,c]) %*% Sinv %*% matZ_ij[,c] )
          }
          for(c in 2:(t_ij-1)){
            S2_jb <- S2_jb + wwij * ( t(matZ_ij[,c]) %*% Sinv %*% matZ_ij[,c] )
          }
          S2_j <- S2_ja + S2_jb
        }
      }
    }
    Fa <- Fa + wi * S1_j
    Fb <- Fb + wi * S2_j
  }
  ### stage 1 estimate of alpha
  alpha0 <- ( Fb - sqrt( ( Fb - 2 * Fa ) * ( Fb + 2 * Fa ) ) ) / ( 2 * Fa )
  return(alpha0)
}