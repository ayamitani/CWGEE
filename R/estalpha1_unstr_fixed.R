#----------------------------------------------------------------------------
# function to estimate stage 1 alpha for unstructured correlation structure
# for fixed number of categories
#----------------------------------------------------------------------------

#' @export

estalpha1_unstr_fixed <- function(mdat, id, Z, nresponse){
  
  Zmat <- 0
  
  for(i in id){
    
    csizei <- nlevels(as.factor(mdat[mdat$cluster == i,]$unit))
    wi <- 1/csizei
    cveci <- unique(mdat[mdat$cluster == i,]$unit)
    
    Zmati <- 0
    
    for (j in cveci){
      
      Z_ij <- Z[mdat$cluster == i & mdat$unit == j]
      Zmati <- Zmati + Z_ij %*% t(Z_ij)
      
    }
    
    Zmat <- Zmat + wi * Zmati
    
  }
  
  delta_k <- delta0 <- diag(rep(1, dim(Zmat)[1]))
  #iter <- 0
  diffR <- 1
  
  while (sum(abs(diffR)) > 0.000001){
    
    delta <- diag( delta_k ^ (1/2) %*% Zmat %*% delta_k ^ (1/2) ) ^ (1/2)
    delta_k1 <- diag( delta )
    diffR <- delta_k1 - delta_k
    delta_k <- delta_k1 
    
    #iter <- iter + 1
    #print(iter)
    
  }
  
  R1 <- solve( delta_k ^ (1/2) ) %*% ( delta_k ^ (1/2) %*% Zmat %*% delta_k ^ (1/2) ) ^ (1/2) %*% solve( delta_k ^ (1/2) )
  
  e <- rep(1, length = dim(R1)[1])
  v_m <- solve( hadamard.prod(R1, R1) ) %*% e
  
  diagv_m <- matrix(0, dim(R1)[1], dim(R1)[1])
  diag(diagv_m) <- v_m
  
  #R <- R1 %*% diagv_m %*% R1
  
  Zmathat <- Zmat/length(id)
  ZZ <- matrix(0, dim(Zmat)[1], dim(Zmat)[1])
  diag(ZZ) <- diag(Zmathat)
  ZZinv <- solve(ZZ)
  R <- (ZZinv) ^ (1/2) %*% Zmathat %*% (ZZinv) ^ (1/2)
  
}


#outc <- estalpha1_unstr_fixed(mdat=out[[1]], id =out[[2]], Z=out[[3]], ncategories = out[[4]])
