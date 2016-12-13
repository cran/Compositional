################################
#### Dirichlet regression for compositional data
#### Tsagris Michail 6/2014
#### mtsagris@yahoo.gr
################################

diri.reg <- function(y, x, plot = TRUE, xnew = NULL) {

  dm <- dim(y)
  n <- dm[1]  ## sample size
  ## the design matrix is created
  x <- model.matrix(y ~ ., data.frame(x) )
  d <- dm[2] - 1  ## dimensionality of the simplex
  z <- log(y)

    dirireg <- function(param, z, x, n, d) {
      phi <- exp( param[1] )  ## this avoids negative values in phi
      para <- param[-1]
      be <- matrix(para, ncol = d)  ## puts the beta parameters in a matrix
      mu1 <- cbind( 1, exp(x %*% be) )
      ma <- mu1 / rowSums(mu1)  ## the fitted values
      ba <- phi * ma
      - n * lgamma(phi) + sum( lgamma(ba) ) - sum( z * (ba - 1 ) )
    }

  runtime <- proc.time()
  rla <- z[, -1] - z[, 1]   ##  log(y[, -1] / y[, 1])  ## additive log-ratio transformation
  ini <- as.vector( coef( lm.fit(x, rla) ) )  ## initial values
  ## based on the logistic normal
  ## the next lines optimize the dirireg function and
  ## estimate the parameter values
  el <- NULL
  options(warn = -1)
  qa <- nlm(dirireg, c(3, ini), z = z, x = x, n = n, d = d)
  el[1] <- -qa$minimum
  qa <- nlm(dirireg, qa$estimate, z = z, x = x, n = n, d = d)
  el[2] <- -qa$minimum
  vim <- 2

  while (el[vim] - el[vim - 1] > 1e-06) {
    ## the tolerance value can of course change
    vim <- vim + 1
    qa <- nlm(dirireg, qa$estimate, z = z, x = x, n = n, d = d)
    el[vim] <- -qa$minimum
  }

  qa <- nlm(dirireg, qa$estimate, z = z, x = x, n = n, d = d, hessian = TRUE)
  log.phi <- qa$estimate[1]
  para <- qa$estimate[-1]  ## estimated parameter values
  be <- matrix(para, ncol = d)  ## matrix of the betas
  colnames(be) <- colnames(y[, -1])  ## names of the betas
  seb <- sqrt( diag( solve(qa$hessian) ) )  ## std of the estimated betas
  std.logphi <- seb[1]  ## std of the estimated log of phi
  seb <- matrix(seb[-1], ncol = d)  ## std of the estimated betas

  if ( !is.null( colnames(y) ) ) {
    colnames(seb) <- colnames(y[, -1])
  } else  colnames(seb) <- paste("Y", 1:d, sep = "")

  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew) )
    mu <- cbind( 1, exp(xnew %*% be) )
    est <- mu / Rfast::rowsums(mu)

  } else {
    mu <- cbind( 1, exp(x %*% be) )
    est <- mu / Rfast::rowsums(mu)  ## fitted values
    lev <- ( exp(log.phi) + 1 ) * Rfast::rowsums( (y - est)^2 / mu )

    if ( plot ) {
      plot(1:n, lev, main = "Influence values", xlab = "Observations",
      ylab = expression( paste("Pearson ", chi^2, "statistic") ) )
      lines(1:n, lev, type = "h")
      abline(h = qchisq(0.95, d), lty = 2, col = 2)
    }

  }

  runtime <- proc.time() - runtime

  rownames(be)  <- colnames(x)
  if  ( !is.null(seb) ) rownames(seb) <- colnames(x)

  list(runtime = runtime, loglik = -qa$minimum, phi = exp(log.phi), log.phi = log.phi,
  std.logphi = std.logphi, be = be, seb = seb, lev = lev, est = est)

}
