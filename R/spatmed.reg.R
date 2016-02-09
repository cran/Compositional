################################
#### Spatial median regression
#### Tsagris Michail 10/2014
#### Biman Chakraborty (2003) On multivariate quantile regression
#### Journal of Statistical Planning and Inference
#### http://www.stat.nus.edu.sg/export/sites/dsap/research/documents/tr01_2000.pdf
#### mtsagris@yahoo.gr
################################

spatmed.reg <- function(y, x, xnew = NULL) {
  ## y contains the dependent variables
  ## x contains the independent variable(s)
  y <- as.matrix(y)
  x <- as.matrix(x)
  d <- ncol(y)  ## dimensionality of y
  x <- cbind(1, x)  ## add the constant term
  p <- ncol(x)  ## dimensionality of x
  z <- list(y = y, x = x)
  ## medi is the function to perform median regression
  medi <- function(beta, z) {
    y <- z$y
    x <- z$x
    p <- ncol(x)
    be <- matrix(beta, nrow = p)
    est <- x %*% be
    sum( sqrt( rowSums((y - est)^2) ) )
  }
  ## we use nlm and optim to obtain the initial beta coefficients
  ini <- matrix(nrow = p, ncol = d)
  for (i in 1:d)  ini[, i] <- coef( quantreg::rq(y[, i] ~ x[, -1]) )
  ini <- as.vector(ini)
  qa <- nlm(medi, ini, z = z, iterlim = 10000)
  qa <- optim(qa$estimate, medi, z = z, control = list(maxit = 20000))
  qa <- optim(qa$par, medi, z = z, control = list(maxit = 20000),
  hessian = TRUE)
  beta <- matrix( qa$par, ncol = ncol(y) )
  if ( is.null(xnew) ) {
    est = x %*% beta
  } else {
    xnew <- cbind(1, xnew)
    xnew <- as.matrix(xnew)
    est <- xnew %*% beta
  }
  sb <- sqrt( diag( solve(qa$hessian) ) )
  sb <- matrix( sb, ncol = ncol(y) )
  if ( is.null(colnames(y)) ) {
    colnames(sb) <- colnames(beta) <- paste("Y", 1:d, sep = "")
  } else  colnames(sb) <- colnames(beta) <- colnames(y)
  if ( is.null(colnames(x)) ) {
    rownames(beta) <- rownames(sb) <- c( "Intercept",
  paste("x", 1:c(p - 1), sep = "") )
  } else  rownames(beta) <- rownames(sb) <- c( "Intercept", colnames(x)[-1] )
  if (d == 1)  est <- as.vector(est)
  list(beta = beta, seb = sb, est = est)
}
