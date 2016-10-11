################################
#### Principal components regression
#### Tsagris Michail 12/2013
#### mtsagris@yahoo.gr
#### References: Jolliffe I.T. (2002)
#### Principal Component Analysis p. 167-188.
################################

pcr <- function(y, x, k = 1, xnew = NULL) {
  ## xnew is the new independent variables values
  ## whose values of y you want to estimate
  ## by default xnew is the x, so you will get the fitted values
  ## y is the univariate dependent variable
  ## x contains the independent variables
  ## k shows the number of components to keep

  x <- as.matrix(x)
  y <- as.vector(y)
  m <- mean(y)
  y <- y - m  ## standardize the dependent variable
  n <- dim(x)[1]
  p <- dim(x)[2]
  mx <- Rfast::colmeans(x) 
  s <- Rfast::colVars(x, std = TRUE) 
  x <- ( t(x) - m )/ s  ## standardise the x values
  x <- t(x)

  eig <- eigen( crossprod(x) )  ## eigen analysis of the design matrix
  values <- eig$values  ## eigenvalues
  per <- cumsum( values / sum(values) )  ## cumulative proportion of each eigenvalue
  vec <- eig$vectors  ## eigenvectors, or principal components
  z <- x %*% vec  ## PCA scores
  mod <- lm.fit(x, y)  ## lm.fit is an internal of lm and thus faster
  sigma <- sum( mod$residuals^2 ) / (n - p - 1)  ## estimated variance
  zzk <- crossprod( z[, 1:k] )
  A <- vec[, 1:k] %*% chol2inv( chol(zzk) )
  b <- A %*% crossprod( z[, 1:k], y )
  ## b is the PCA based coefficients

  mse <- r2 <- NULL
  if ( !is.null(xnew) ) {
    xnew <- matrix(xnew, ncol = p)
    nu <- dim(xnew)[1]
    xnew <- ( t(xnew) - mx) / s  ## standardize the xnew values
    xnew <- t(xnew)
    est <- as.vector( m + xnew %*% b )  ## predicted values for PCA model
  } else {
    est <- as.vector( m + x %*% b )  ## fitted values for PCA model
    mse <- sum( (y + m - est)^2 ) / (n - k)  ## mean squared error of PCA model
    r2 <- 1 - (n - 1)/(n - k - 1) * ( 1 - cor(y + m, est)^2 )
  }

  ## rs is the adjusted R squared for the PCA model
  va <- sigma * tcrossprod( A, vec[, 1:k] )
  ## va is the covariance matrix of the parameters
  ## of the parameters of the PCA model
  vara <- sqrt( diag(va) )  ## standard errors of coefficients of PCA model

  param <- cbind(b, vara)
  colnames(param) <- c("beta", "std.error")
  list(parameters = param, mse = mse, adj.rsq = r2, per = per[k], est = est)

}
