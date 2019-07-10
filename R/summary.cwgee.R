#------------------------------------
# Summary
#------------------------------------

#' @export

summary.cwgee <- function(object, ...) {
  robust.se <- sqrt(diag(object$robust.variance))
  robust.z <- coef(object)/robust.se
  pvalue <- 2 * (1 - pnorm(abs(robust.z)))
  TAB <- cbind(Estimate = coef(object), 
               Robust_S.E. = robust.se, 
               Robust_z = robust.z, 
               p.value = pvalue)
  TAB <- round(TAB, 5)
  colnames(TAB) <- c("Estimate", "san.se", "san.z", "Pr(>|san.z|)")
  rownames(TAB) <- object$coef.names
  res <- list(coefficients = TAB, 
              call = object$call,
              corr.str = object$corr.str,
              corr.matrix = object$corr.matrix,
              corr.coef = object$alpha,
              time.str = object$time.str,
              niter = object$niter)
  class(res) <- "summary.cwgee"
  res
}