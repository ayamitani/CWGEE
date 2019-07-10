#------------------------------------------------------------
# function to estimate stage 1 alpha for exch structure
#------------------------------------------------------------


estalpha1_exch <- function(mdat, id, Z, Sinv, ncategories1){
  
  match.call()
  
  alphafun <- function(alpha){
    
    GG1 <- GG2 <- 0
    for (i in id){
      csizei <- nlevels(as.factor(mdat[mdat$id == i,]$cluster.var))
      wi <- 1/csizei
      cveci <- unique(mdat[mdat$id == i,]$cluster.var)
      GG1j <- GG2j <- 0
      for (j in cveci){
        t_ij <- nlevels(as.factor(mdat[mdat$id == i & mdat$cluster.var == j,]$time.var))
        wwij <- 1 / t_ij
        Z_ij <- Z[mdat$id == i & mdat$cluster.var == j]
        matZ_ij <- matrix(Z_ij, nrow = ncategories1) 
        if(t_ij > 1){
          g1 <- vector(length = t_ij)
          for(t in 1:t_ij) g1[t] <- t(matZ_ij[,t]) %*% Sinv %*% matZ_ij[,t]
          G1 <- sum(g1)
          
          g2 <- vector()
          for(t in 1:(t_ij - 1)){
            for(tt in (t+1):t_ij){
              g2 <- c(g2, t(matZ_ij[,t]) %*% Sinv %*% matZ_ij[,tt])
            }
          }
          G2 <- sum(g2)
          
          denom <- ( 1 + ( t_ij - 1 ) * alpha ) ^ 2
          num1 <- alpha ^ 2 * ( t_ij - 1 ) * ( t_ij - 2 ) + 2 * alpha * ( t_ij - 1 )
          num2 <- ( 1 + alpha ^ 2 * ( t_ij - 1 ) )
        }
        GG1j <- GG1j + wwij * ( G1 * num1 ) / denom
        GG2j <- GG2j + wwij * ( G2 * num2 ) / denom
      }
      GG1 <- GG1 + wi * GG1j
      GG2 <- GG2 + wi * GG2j
    }
    GG1 - 2 * GG2
  }
  
  ### stage 1 estimate of alpha
  alpha0 <- uniroot(alphafun, c(0,1), tol = 1e-10, extendInt = "yes")$root
  return(alpha0)
}