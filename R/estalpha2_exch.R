#------------------------------------------------------------
# function to estimate stage 2 alpha for exch structure
#------------------------------------------------------------

estalpha2_exch <- function(alpha0, mdat, id){
  
  match.call()
  
  alphapart1 <- alphapart2 <- 0
  for (i in id){
    csizei <- nlevels(as.factor(mdat[mdat$id == i,]$cluster.var))
    wi <- 1/csizei
    cveci <- unique(mdat[mdat$id == i,]$cluster.var)
    alphapart1j <- alphapart2j <- 0
    for (j in cveci){
      t_ij <- nlevels(as.factor(mdat[mdat$id == i & mdat$cluster.var == j,]$time.var))
      if(t_ij > 1){
        alphapart1num <- alpha0 * ( t_ij - 1 )* ( alpha0 * (t_ij - 2) + 2 )
        alphapart2num <- ( t_ij - 1 ) * ( 1 + alpha0 ^ 2 * (t_ij - 1) )
        alphaden <- ( 1 + alpha0 * ( t_ij - 1 ) ) ^ 2
        
        alphapart1j <- alphapart1j + alphapart1num / alphaden
        alphapart2j <- alphapart2j + alphapart2num / alphaden
      }
    }
    alphapart1 <- alphapart1 + wi * alphapart1j
    alphapart2 <- alphapart2 + wi * alphapart2j
  }
  alpha <- alphapart1 / alphapart2
  return(alpha)
}