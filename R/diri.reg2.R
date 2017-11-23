################################
#### Dirichlet regression with covariates on phi
####  2/2015
#### mtsagris@yahoo.gr
################################
diri.reg2 <- function(y, x, xnew = NULL) {

  n <- dim(y)[1]  ## sample size
  x <- model.matrix(y ~ ., data.frame(x) )
  p <- dim(x)[2]  ## dimensionality of x
  d <- dim(y)[2] - 1
  ly <- Rfast::Log(y)  ## dimensionality of the simplex

    dirireg2 <- function(param) {
      ## param contains the parameter values
      phipar <- param[1:p]
      para <- param[ -c(1:p) ]
      phi <- exp( x %*% phipar )  ## phi is a function of the covariates
      be <- matrix(para, nrow = p)  ## puts the beta parameters in a matrix
      mu1 <- cbind( 1, exp(x %*% be) )
      ma <- mu1 / rowSums(mu1)  ## the fitted values
      ba <- as.vector(phi) * ma
      - sum( lgamma(phi) ) + sum( lgamma(ba) ) - sum( ly * (ba - 1) )
    }

  runtime <- proc.time()
  rla <- ly[, -1] - ly[, 1]    ##  log( y[, -1] / y[, 1] )  ## additive log-ratio transformation
  ini <- as.vector( lm.fit(x, rla)$coefficients )  ## initial values
  ## based on the logistic normal
  ## the next lines optimize the dirireg2 function and
  ## estimate the parameter values

  el <- NULL
  qa <- nlm(dirireg2, c(rnorm(p, 0, 0.1), as.vector( t(ini) ) ) )
  el1 <-  -qa$minimum
  qa <- nlm(dirireg2, qa$estimate)
  el2 <- -qa$minimum
  while (el2 - el1 > 1e-06) {
    ## the tolerance value can of course change
    el1 < -el2
    qa <- nlm(dirireg2, qa$estimate)
    el2 <-  -qa$minimum
  }

  qa <- optim(qa$estimate, dirireg2, hessian = TRUE)
  phipar <- qa$par[1:p]
  be <- matrix(qa$par[-c(1:p)], nrow = p)  ## matrix of the betas
  phi <- as.numeric( exp(x %*% phipar) )  ## estimated beta parameters of phi
  s <- sqrt( diag( solve(qa$hessian) ) )  ## std of the estimated parameters
  std.phi <- s[1:p]  ## std of the estimated beta parameters of the phi
  seb <- matrix( s[-c(1:p)], ncol = d )  ## std of the estimated betas
  V <- solve(qa$hessian)  ## covariance matrix of the parameters
  runtime <- proc.time() - runtime

  if ( !is.null( colnames(y) ) ) {
    colnames(be) <- colnames(seb) <- colnames(y[, -1])
  } else  colnames(beta) <- colnames(seb) <- paste("Y", 1:d, sep = "")
  
  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew) )
    mu <- cbind( 1, exp(xnew %*% be) )
    est <- mu / Rfast::rowsums(mu)
  } else  est <- NULL

  rownames(be)  <- colnames(x)
  if  ( !is.null(seb) ) rownames(seb) <- colnames(x)
  list(runtime = runtime, loglik = -qa$value, phipar = phipar,
       std.phi = std.phi, be = be, seb = seb, sigma = V, phi = phi, est = est)
}
