#------------------------------------------------------------
# function for iterative CWGEE
#------------------------------------------------------------

#' @export

ordCWGEEgen <- function(mdat, beta0, nbeta, id, Y_var, ncategories1, ncategories, X_mat, time.str, Z){
  
  match.call()
  beta <- beta0
  bdiff <- 1
  iter <- 0
  while(sum(abs(bdiff))>.00000001){
    DWZ <- matrix(0, nrow = nbeta)
    DWD <- matrix(0, ncol = nbeta, nrow = nbeta)
    DWZZWD <- matrix(0, ncol = nbeta, nrow = nbeta)
    
    ### matrix of correlation between categories
    S_mat <- smat.work(beta[1:ncategories1], ncategories)
    S <- S_mat$S_mat
    Sinv <- solve(S)
    
    if (time.str == "ind")
      alpha0 <- 0
    if (time.str == "ar1")
      alpha0 <- estalpha1_ar1(mdat = mdat, id = id, Z = Z, Sinv = Sinv, ncategories1 = ncategories1)
    if (time.str == "exch")
      alpha0 <- estalpha1_exch(mdat = mdat, id = id, Z = Z, Sinv = Sinv, ncategories1 = ncategories1)
    
    for (i in id){
      csizei <- nlevels(as.factor(mdat[mdat$id == i,]$cluster.var))
      wi <- 1/csizei
      cveci <- unique(mdat[mdat$id == i,]$cluster.var)
      DWZj <- matrix(0, nrow = nbeta)
      DWDj <- matrix(0, ncol = nbeta, nrow = nbeta)
      for(j in cveci){
        t_ij <- nlevels(as.factor(mdat[mdat$id == i & mdat$cluster.var == j,]$time.var))
        wwij <- 1 / t_ij
        C <- corr.str(d = t_ij, time.str = time.str, rho = alpha0)
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
    Xb <- X_mat %*% beta
    U <- exp(Xb) / (1 + exp(Xb))
    var <- U * (1 - U)
    Z <- (Y_var - U) / sqrt(var)
    iter <- iter + 1
    #print(iter)
  }
  
  fit <- list()
  fit$coefficients <- beta
  fit$robust.variance <- covbeta
  fit$robust.se <- sqrt(diag(covbeta))
  fit$wald.chisq <- (beta/sqrt(diag(covbeta)))^2
  fit$p.value <- 1 - pchisq(fit$wald.chisq, df = 1)
  fit$alpha <- alpha0
  fit$niter <- iter
  fit
  
}