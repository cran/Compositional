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
  m <- mean(y)
  y <- y - m  ## standardize the dependent variable
  p <- dim(x)[2]
  if (k > p)   k <- p
  eig <- prcomp(x, center = TRUE, scale = TRUE)
  values <- eig$sdev^2
  per <- cumsum( values / sum(values) )  ## cumulative proportion of each eigenvalue
  vec <- eig$rotation[, 1:k, drop = FALSE]
  z <- x %*% vec  ## PCA scores
  zk <- Rfast::colsums(z^2)
  zzk <- matrix(0, k, k)
  diag(zzk) <- zk
  A <- vec %*% chol2inv( chol(zzk) )
  be <- A %*% crossprod( z, y )  ## PCA based coefficients
  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- matrix(xnew, ncol = p)
    est <- as.vector( m + xnew %*% be )  ## predicted values for PCA model
  }
  rownames(be) <- colnames(x)
  list(parameters = be, per = per[k], est = est)
}
