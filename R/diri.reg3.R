diri.reg3 <- function(y, x, xnew = NULL) {

  dm <- dim(y)
  n <- dm[1]  ## sample size
  ## the design matrix is created
  x <- model.matrix(y ~ ., data.frame(x) )
  d <- dm[2]  ## dimensionality of the simplex
  z <- Rfast::Log(y)

  dirireg3 <- function(param, z, x, d) {
    be <- matrix(param, ncol = d)  ## puts the beta parameters in a matrix
    a <- exp(x %*% be)
    phi <- Rfast::rowsums(a)
    - sum( lgamma(phi) ) + sum( lgamma(a) ) - sum( z * (a - 1 ) )
  }

  runtime <- proc.time()
  ini.beta <- NULL
  for (i in 1:d)  ini.beta <- c(ini.beta, Rfast2::gammareg(y[, 1], x[, -1])$be)
  ## based on the logistic normal
  ## the next lines optimize the dirireg function and
  ## estimate the parameter values
  suppressWarnings({
    qa <- optim( ini.beta, dirireg3, z = z, x = x, d = d, method = "BFGS", control = list(maxit = 5000) )
    qa <- optim( qa$par, dirireg3, z = z, x = x, d = d, method = "BFGS", control = list(maxit = 5000)  )
    qa <- optim(qa$par, dirireg3, z = z, x = x, d = d, method = "BFGS", control = list(maxit = 5000), hessian = TRUE)
  })
  be <- matrix(qa$par, ncol = d)  ## matrix of the betas
  colnames(be) <- colnames(y)  ## names of the betas
  seb <- sqrt( diag( solve(qa$hessian) ) )  ## std of the estimated betas
  seb <- matrix(seb, ncol = d)  ## std of the estimated betas

  if ( !is.null( colnames(y) ) ) {
    colnames(seb) <- colnames(y)
  } else  colnames(seb) <- paste("Y", 1:d, sep = "")

  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew) )
    mu <- exp(xnew %*% be)
    est <- mu / Rfast::rowsums(mu)
    lev <- NULL
  } else {
    mu <- exp(x %*% be)
    est <- mu / Rfast::rowsums(mu)  ## fitted values
  }

  runtime <- proc.time() - runtime
  rownames(be)  <- colnames(x)
  if  ( !is.null(seb) ) rownames(seb) <- colnames(x)
  list(runtime = runtime, loglik = -qa$value, be = be, seb = seb, est = est)
}
