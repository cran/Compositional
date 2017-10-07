################################
#### Multivariate (and univariate) ridge regression
#### Tsagris Michail 10/2014
#### mtsagris@yahoo.gr
################################
ridge.reg <- function(y, x, lambda, B = 1, xnew = NULL) {
  ## xnew is the new independent variables values
  ## whose values of y you want to estimate
  ## by default xnew is the x, so you will get the fitted values
  ## y is a real valued vector
  ## x contains the independent variable(s)
  ## lambda is the ridge regularization parameter
  ## if lambda=0, the classical multivariate regression is implemented
  ## B is for bootstrap estimation of the standard errors of the betas
  ## if pred is TRUE it means that you want to predict new y
  ## but if xnew is x (by default), the pred is not important
  ## the pred is important if xnew is not x
  if ( all( y > 0 & y < 1 ) )  y <- log(y / ( 1 - y) ) ## logistic normal
  n <- length(y)  ## sample size
  p <- dim(x)[2]  ## dimensionality of x
  my <- sum(y) / n
  yy <- y - my   ## center the dependent variables
  mx <- Rfast::colmeans(x)
  s <- Rfast::colVars(x, std = TRUE)
  xx <- t( ( t(x) - mx )/s ) ## standardize the independent variables
  lamip <- lambda * diag(p)
  xtx <- crossprod(xx)
  W <- solve( xtx + lamip )
  beta <- W %*% crossprod(xx, yy)
  est <- as.vector( xx %*% beta + my )
  va <- Rfast::Var(y - est) * (n - 1) / (n - p - 1)
  # vab <- kronecker(va, W %*% xtx %*% W  )
  # seb <- matrix( sqrt( diag(vab) ), nrow = p )
  vab <- va * mahalanobis(W, numeric(p), xtx, inverted = TRUE)
  seb <- sqrt( vab )
  if (B > 1) { ## bootstrap estimation of the standard errors
    be <- matrix(nrow = B, ncol = p )
    for ( i in 1:B) {
      id <- sample(1:n, n, replace = TRUE)
      yb <- yy[id, ]     ;     xb <- xx[id, ]
      be[i, ] <- solve( crossprod(xb) + lamip, crossprod(xb, yb) )
    }
    seb <- Rfast::colVars(be, std = TRUE) ## bootstrap standard errors of betas
  }
  ## seb contains the standard errors of the coefficients
  if ( is.null( colnames(x) ) ) {
    names(seb) <- paste("X", 1:p, sep = "")
    names(beta) <- paste("X", 1:p, sep = "")
  } else  names(seb) <- names(beta) <- colnames(x)

  if ( !is.null(xnew) ) {
    xnew <- matrix(xnew, ncol = p)
    xnew <- t( ( t(xnew) - mx ) / s ) ## scale the xnew values
    est <- as.vector( xnew %*% beta + my )
  } else est <- est

  list(beta = beta, seb = seb, est = est)
}
