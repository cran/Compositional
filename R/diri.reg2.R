################################
#### Dirichlet regression with covariates on phi
####  2/2015
#### mtsagris@yahoo.gr
################################

diri.reg2 <- function(y, x, xnew = NULL) {
  ## y is the compositional data
  y <- as.matrix(y)
  y <- y / as.vector( Rfast::rowsums(y) )
  n <- nrow(y)  ## sample size
  mat <- model.matrix(y ~ ., as.data.frame(x) )
  x <- as.matrix(mat[1:n, ])  ## the desing matrix is created
  p <- ncol(x)  ## dimensionality of x
  ## the line above makes sure y is compositional data and
  ## then the unit vector is added to the design matrix
  d <- ncol(y) - 1
  z <- list(y = log(y), x = x)  ## dimensionality of the simplex

    dirireg2 <- function(param, z = z) {
      ## param contains the parameter values
      ## z contains the compositional data and independent variable(s)
      y <- z$y
      x <- z$x
      ## y is the compositional data and x the independent variable(s)
      p <- ncol(x)  ## dimensionality of x
      n <- nrow(y)  ## sample size
      d <- ncol(y) - 1  ## dimensionality of  the simplex
      phipar <- param[1:p]
      para <- param[ -c(1:p) ]
      phi <- exp( x %*% phipar )  ## phi is a function of the covariates
      be <- matrix(para, ncol = d)  ## puts the beta parameters in a matrix
      mu1 <- cbind( 1, exp(x %*% be) )
      ma <- mu1 / rowSums(mu1)  ## the fitted values
      ba <- as.vector(phi) * ma
      l <-  -( sum( lgamma(phi) ) - sum( lgamma(ba) ) +   sum( y * (ba - 1) ) )
      ## l is the log-likelihood
      l
    }

  runtime <- proc.time()
  rla <- log( y[, -1] / y[, 1] )  ## additive log-ratio transformation
  ini <- as.vector( coef( lm.fit(x, rla) ) )  ## initial values
  ## based on the logistic normal
  ## the next lines optimize the dirireg2 function and
  ## estimate the parameter values

  el <- NULL
  qa <- nlm(dirireg2, c(rnorm(p, 0, 0.1), as.vector( t(ini) ) ), z = z)
  el[1] <- -qa$minimum
  qa <- nlm(dirireg2, qa$estimate, z = z)
  el[2] <- -qa$minimum
  vim <- 2
  while (el[vim] - el[vim - 1] > 1e-04) {
    ## the tolerance value can of course change
    vim <- vim + 1
    qa <- nlm(dirireg2, qa$estimate, z = z)
    el[vim] <- -qa$minimum
  }

  qa <- nlm(dirireg2, qa$estimate, z = z, hessian = TRUE)
  phipar <- qa$estimate[1:p]
  para <- qa$estimate[-c(1:p)]  ## estimated parameter values
  beta <- matrix(para, ncol = d)  ## matrix of the betas
  mu1 <- cbind( 1, exp(x %*% beta) )
  ma <- mu1 / rowSums(mu1)  ## fitted values
  phi <- as.numeric( exp(x %*% phipar) )  ## estimated beta parameters of phi
  s <- sqrt( diag( solve(qa$hessian) ) )  ## std of the estimated parameters
  std.phi <- s[1:p]  ## std of the estimated beta parameters of the phi
  seb <- matrix( s[-c(1:p)], ncol = d )  ## std of the estimated betas
  V <- solve(qa$hessian)  ## covariance matrix of the parameters

  runtime <- proc.time() - runtime

  if ( !is.null( colnames(y) ) ) {
    colnames(beta) <- colnames(seb) <- colnames(y[, -1])
  } else  colnames(beta) <- colnames(seb) <- paste("Y", 1:d, sep = "")

  if ( !is.null(xnew) ) {
    xnew <- cbind(1, xnew)
    xnew <- as.matrix(xnew)
    mu <- cbind( 1, exp(xnew %*% beta) )
    est <- mu / as.vector( Rfast::rowsums(mu) )

  } else {
    mu <- cbind( 1, exp(x %*% beta) )
    est <- mu / as.vector( Rfast::rowsums(mu) ) ## fitted values
  }


  if ( is.null(colnames(x)) ) {
    p <- ncol(x) - 1
    rownames(beta) <- c("constant", paste("X", 1:p, sep = "") )
    if ( !is.null(seb) )  rownames(seb) <- c("constant", paste("X", 1:p, sep = "") )
  } else {
    rownames(beta)  <- c("constant", colnames(x)[-1] )
    if  ( !is.null(seb) ) rownames(seb) <- c("constant", colnames(x)[-1] )
  }

  list(runtime = runtime, loglik = -qa$minimum, phipar = phipar,
       std.phi = std.phi, beta = beta, seb = seb, sigma = V, phi = phi, est = est)

}
