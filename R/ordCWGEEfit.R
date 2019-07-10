#------------------------------------------------------------
# function to fit CWGEE with alpha
#------------------------------------------------------------

#' @export

ordCWGEEfit <- function(mdat, beta0, nbeta, id, Y_var, ncategories1, ncategories, X_mat, time.str, Z, alpha){
  
  match.call()
  beta <- beta0
  DWZ <- matrix(0, nrow = nbeta)
  DWD <- matrix(0, ncol = nbeta, nrow = nbeta)
  DWZZWD <- matrix(0, ncol = nbeta, nrow = nbeta)
  
  ### matrix of correlation between categories
  S_mat <- smat.work(beta[1:ncategories1], ncategories)
  S <- S_mat$S_mat
  Sinv <- solve(S)
  
  for (i in id){
    csizei <- nlevels(as.factor(mdat[mdat$id == i,]$cluster.var))
    wi <- 1/csizei
    cveci <- unique(mdat[mdat$id == i,]$cluster.var)
    DWZj <- matrix(0, nrow = nbeta)
    DWDj <- matrix(0, ncol = nbeta, nrow = nbeta)
    for(j in cveci){
      t_ij <- nlevels(as.factor(mdat[mdat$id == i & mdat$cluster.var == j,]$time.var))
      wwij <- 1 / t_ij
      C <- corr.str(d = t_ij, time.str = time.str, rho = alpha)
      R <- C %x% S
      y <- as.matrix( Y_var[mdat$id == i & mdat$cluster.var == j] )
      x <- as.matrix( X_mat[mdat$id == i & mdat$cluster.var == j,] )
      xb <- x %*% beta
      u <- exp( xb ) / ( 1 + exp( xb ) ) 
      dudb <- exp( xb ) / ( 1 + exp( xb ) ) ^ 2
      D <- x[,1] * dudb
      for (p in 2:ncol(x)){
        D <- cbind(D, x[,p] * dudb)
      }
      vv <- u * ( 1 - u )
      V <- matrix(ncol = length(vv), nrow = length(vv), 0)
      diag(V) <- vv
      W <- V ^ (1/2) %*% R %*% V ^ (1/2)
      invW <- solve(W)
      DWZj <- DWZj + wwij * t(D) %*% invW %*% ( y - u ) 
      DWDj <- DWDj + wwij * t(D) %*% invW %*% D
    }
    DWZ <- DWZ + wi * DWZj
    DWD <- DWD + wi * DWDj
    DWZZWD <- DWZZWD + ( wi * DWZj ) %*% t( wi * DWZj )
  }
  
  invDWD <- solve(DWD)
  bdiff <- invDWD %*% DWZ
  beta <- beta + bdiff
  covbeta <- invDWD %*% DWZZWD %*% invDWD
  
  out <- list()
  out$call <- match.call()
  out$coefficients <- beta
  out$robust.variance <- covbeta
  out$robust.se <- sqrt(diag(covbeta))
  out$wald.chisq <- (beta/sqrt(diag(covbeta)))^2
  out$p.value <- 1 - pchisq(out$wald.chisq, df = 1)
  out$alpha <- alpha
  class(out) <- "cwgee"
  out
  
}