################################
#### Multivariate (and univariate) ridge regression
#### Tsagris Michail 10/2014
#### mtsagris@yahoo.gr
################################
ridge.reg <- function(y, x, lambda, B = 1, xnew = NULL) {
  ## xnew is the new independent variables values
  ## whose values of y you want to estimate
  ## by default xnew is NULL
  ## y is a real valued vector
  ## x contains the independent variable(s)
  ## lambda is the ridge regularization parameter
  ## if lambda=0, the classical regression is implemented
  ## B is for bootstrap estimation of the standard errors of the betas
  if ( all( y > 0  &  y < 1 ) )  y <- log(y / ( 1 - y) ) ## logistic normal
  n <- length(y)  ## sample size
  p <- dim(x)[2]  ## dimensionality of x
  my <- sum(y) / n
  yy <- y - my   ## center the dependent variables
  lamip <- lambda * diag(p)
  xtx <- crossprod(x)
  W <- solve( xtx + lamip )
  beta <- W %*% crossprod(x, yy)
  est <- as.vector( x %*% beta + my )
  va <- Rfast::Var(y - est) * (n - 1) / (n - p - 1)
  vab <- va * mahalanobis(W, numeric(p), xtx, inverted = TRUE)
  seb <- sqrt( vab )
  
  if (B > 1) { ## bootstrap estimation of the standard errors
    be <- matrix(nrow = B, ncol = p )
    for ( i in 1:B) {
      id <- Rfast2::Sample.int(n, n, replace = TRUE)
      yb <- yy[id, ]     ;     xb <- x[id, ]
      be[i, ] <- solve( crossprod(xb) + lamip, crossprod(xb, yb) )
    }
    seb <- Rfast::colVars(be, std = TRUE) ## bootstrap standard errors of betas
  }
  
  ## seb contains the standard errors of the coefficients
  if ( is.null( colnames(x) ) ) {
    names(seb) <- paste("X", 1:p, sep = "")
    names(beta) <- paste("X", 1:p, sep = "")
  } else  names(seb) <- names(beta) <- colnames(x)

  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- matrix(xnew, ncol = p)
    est <- as.vector( xnew %*% beta + my )
  }

  list(beta = beta, seb = seb, est = est)
}
