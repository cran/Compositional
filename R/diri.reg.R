diri.reg <- function(y, x, plot = FALSE, xnew = NULL) {

  dm <- dim(y)
  n <- dm[1]  ## sample size
  ## the design matrix is created
  x <- model.matrix(y ~ ., data.frame(x) )
  d <- dm[2] - 1  ## dimensionality of the simplex
  z <- Rfast::Log(y)

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
  ini.beta <- as.vector( Compositional::kl.compreg(y, x[, -1], con = TRUE)$be )
  ini.phi <- sum( Compositional::diri.nr(y)$param )
  ## based on the logistic normal
  ## the next lines optimize the dirireg function and
  ## estimate the parameter values
  el <- NULL
  suppressWarnings({
    qa <- optim( c(ini.phi, ini.beta), dirireg, z = z, x = x, n = n, d = d, control = list(maxit = 5000) )
    qa <- optim( qa$par, dirireg, z = z, x = x, n = n, d = d, control = list(maxit = 5000)  )
    qa <- optim(qa$par, dirireg, z = z, x = x, n = n, d = d, control = list(maxit = 5000), hessian = TRUE)
  })
  log.phi <- qa$par[1]
  be <- matrix(qa$par[-1], ncol = d)  ## matrix of the betas
  colnames(be) <- colnames(y[, -1])  ## names of the betas

  seb <- try( solve( qa$hessian ), silent = TRUE )
  if ( !identical( class(seb), "try-error") ) {
    std.logphi <- seb[1]  ## std of the estimated log of phi
    seb <- matrix( sqrt( diag(seb)[-1] ), ncol = d)  ## std of the estimated betas
    if ( !is.null( colnames(y) ) ) {
      colnames(seb) <- colnames(y[, -1])
    } else  colnames(seb) <- paste("Y", 1:d, sep = "")
    rownames(seb) <- colnames(x)
  } else  seb <- NULL

  est <- NULL
  lev <- NULL

  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew) )
    mu <- cbind( 1, exp(xnew %*% be) )
    est <- mu / Rfast::rowsums(mu)
  } else {
    if ( plot ) {
    mu <- cbind( 1, exp(x %*% be) )
    est <- mu / Rfast::rowsums(mu)  ## fitted values
    lev <- ( exp(log.phi) + 1 ) * Rfast::rowsums( (y - est)^2 / mu )
      plot( 1:n, lev, main = "Influence values", xlab = "Observations",
            ylab = expression( paste("Pearson ", chi^2, "statistic") ),
	        cex.lab = 1.2, cex.axis = 1.2 )
      lines(1:n, lev, type = "h")
      abline(h = qchisq(0.95, d), lty = 2, col = 2)
    }
  }

  runtime <- proc.time() - runtime
  rownames(be)  <- colnames(x)
  if  ( !is.null(seb) ) rownames(seb) <- colnames(x)
  list(runtime = runtime, loglik = -qa$value, phi = exp(log.phi), log.phi = log.phi,
  std.logphi = std.logphi, be = be, seb = seb, lev = lev, est = est)
}
