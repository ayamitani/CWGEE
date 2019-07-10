#' Cluster weighted GEE for ordinal clustered longitudinal data with informative cluster size.
#'
#' Solves the cluster-weighted generalized estimating equations for correlated ordinal
#' responses in clustered longitudinal data assuming a
#' cumulative link logit model for the marginal probabilities
#' using the method of quasi-least squares.
#'
#' The \code{data} must be provided in case level or equivalently in `long' format.
#'
#' @param formula a formula expression as for other regression models.
#' @param data an optional data frame containing the variables provided in
#' \code{formula}, \code{id}, \code{cluster.var} and \code{time.var}.
#' @param id a vector that identifies the clusters.
#' @param cluster.var a vector that identifies the unit within a cluster.
#' @param time.var a vector that identifies the repeated observation of a unit.
#' @param time.str a character string that indicates the temporal working correlation structure.
#' Options include \code{"ind"} for independence, \code{"ar1"} for AR1, and \code{"exch"} for exchangeable.
#'
#' @return Returns an object of the class \code{"cwgee"}. This has components:
#' \item{call}{the matched call.}
#' \item{coefficients}{the estimated regression parameter vector of the marginal model.}
#' \item{coef.names}{the variable name of the coefficients.}
#' \item{robust.variance}{the estimated "robust" covariance matrix.}
#' \item{robust.se}{the estimated "robust" standard errors.}
#' \item{wald.chisq}{the Wald Chi-square test statistic for coefficient estimates.}
#' \item{p.value}{the p-value based on a Wald Chi-square test statistic that no covariates are statistically significant.}
#' \item{alpha}{the estimated temporal correlation coefficient.}
#' \item{niter}{the number of iterations the model took to converge.}
#' \item{time.str}{the temporal working correlation structure assumed for the model.}
#' @author Aya Mitani
#' @examples
#' data(perio)
#' fitmod <- ordCWGEE(formula = cal ~ mets + edu + age + smoking, data = perio,
#' id = subject, cluster.var = tooth, time.var = visit, time.str = "ind")
#' summary(fitmod)
#' @export



ordCWGEE <- function(formula, data, id, cluster.var, time.var, time.str){

  call <- match.call()
  mcall <- match.call(expand.dots = FALSE)
  mf <- match(c("formula", "data", "id", "cluster.var", "time.var"), names(mcall), 0L)
  m <- mcall[c(1L, mf)]
  if (is.null(m$id))
    m$id <- as.name("id")
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  Y <- as.numeric(factor(model.response(m)))
  if (is.null(Y)){
    stop("response variable not found")
  }
  ncategories <- nlevels(factor(Y))
  ncategories1 <- ncategories - 1
  if (ncategories <= 2)
    stop("response variable must have more than 2 categories")
  id <- model.extract(m, "id")
  cluster.var <- model.extract(m, "cluster.var")
  time.var <- model.extract(m, "time.var")
  if (is.null(id)){
    stop("'id' variable not found")
  }
  if (is.null(cluster.var)){
    stop("'cluster.var' variable not found")
  }
  if (is.null(time.var)){
    stop("'time.var' variable not found")
  }
  if (length(id) != length(Y)){
    stop("response variable and 'id' variable are not of same length")
  }
  if (length(cluster.var) != length(Y)){
    stop("response variable and 'cluster.var' variable are not of same length")
  }
  if (length(time.var) != length(Y)){
    stop("response variable and 'time.var' variable are not of same length")
  }
  time.strs <- c("ind", "exch", "ar1")
  strcheck <- match(time.str, time.strs, -1)
  if (strcheck < 1){
    stop("unknown working correlation structure")
  }
  mterms <- attr(m, "terms")
  xvars <- as.character(attr(mterms, "variables"))[-c(1:2)]
  nbeta <- length(xvars) + ncategories1
  names(m) <- c("y", xvars, "id", "cluster.var", "time.var")
  exdata <- ord.expand(scores="y", data=m, subjects="id", categories=ncategories)
  mdat <- as.data.frame(exdata$exdata)
  # model matrix
  mformula <- as.formula(paste("y", "~", paste(c("factor(cuts)-1", paste(xvars, collapse="+")), collapse="+")))
  X_mat <- model.matrix(mformula, data = exdata$exdata)
  mframe <- model.frame(mformula, data = exdata$exdata)
  Y_var <- model.response(mframe)
  unique.id <- unique(id)

  # use GLM to get initial beta estimates
  beta0 <- as.vector(coef(glm(mformula, data = exdata$exdata, family = "binomial")))
  Z <- glm(mformula, data = exdata$exdata, family = "binomial")$residuals
  #return(list(mdat, beta0, nbeta, unique.id, Y_var, ncategories1, ncategories, X_mat, time.str, Z))
  # compute beta and stage 1 alpha estimate
  cwgeeout1 <- ordCWGEEgen(mdat = mdat, beta0 = beta0, nbeta = nbeta, id = unique.id, Y_var = Y_var, ncategories1 = ncategories1, ncategories = ncategories, X_mat = X_mat, time.str = time.str, Z = Z)
  niter <- cwgeeout1$niter
  stage1alpha <- cwgeeout1$alpha
  stage1coefs <- cwgeeout1$coefficients
  # compute stage 2 alpha estimate and final beta estimate
  if (time.str == "ind")
    cwgeeout <- cwgeeout1
  if (time.str == "ar1"){
    alphastage2 <- estalpha2_ar1(alpha0 = stage1alpha)
    cwgeeout <- ordCWGEEfit(mdat = mdat, beta0 = stage1coefs, nbeta = nbeta, id = unique.id, Y_var = Y_var, ncategories1 = ncategories1, ncategories = ncategories, X_mat = X_mat, time.str = time.str, Z = Z, alpha = alphastage2)
  }
  if (time.str == "exch"){
    alphastage2 <- estalpha2_exch(alpha0 = stage1alpha, mdat = mdat, id = id)
    cwgeeout <- ordCWGEEfit(mdat = mdat, beta0 = stage1coefs, nbeta = nbeta, id = unique.id, Y_var = Y_var, ncategories1 = ncategories1, ncategories = ncategories, X_mat = X_mat, time.str = time.str, Z = Z, alpha = alphastage2)
  }

  intnum <- 1:ncategories1
  coef.names <- c(paste0("int", intnum), xvars)

  result <- list()
  result$call <- match.call()
  result$coefficients <- cwgeeout$coefficients
  result$coef.names <- coef.names
  result$robust.variance <- cwgeeout$robust.variance
  result$robust.se <- cwgeeout$robust.se
  result$wald.chisq <- cwgeeout$wald.chisq
  result$p.value <- cwgeeout$p.value
  result$alpha <- cwgeeout$alpha
  result$niter <- niter
  result$time.str <- time.str
  class(result) <- "cwgee"
  result

}
