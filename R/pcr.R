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
  my <- mean(y)
  y <- y - my  ## standardize the dependent variable
  dm <- dim(x)
  n <- dm[1]
  p <- dm[2]
  if (k > p)   k <- p
  m <- Rfast::colmeans(x)
  s <- Rfast::colVars(x, suma = n * m, std = TRUE)
  x <- t(x) - m
  x <- x/s
  x <- t(x)
  eig <- prcomp(x, center = FALSE, scale = FALSE)
  values <- eig$sdev^2
  per <- cumsum( values / sum(values) )  ## cumulative proportion of each eigenvalue
  vec <- eig$rotation[, 1:k, drop = FALSE]
  z <- x %*% vec  ## PCA scores
  zzk <- crossprod(z)
  be <- vec %*% Rfast::spdinv(zzk) %*% crossprod( z, y )  ## PCA based coefficients
  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- matrix(xnew, ncol = p)
    xnew <- t(xnew) - m
    xnew <- t( xnew / s )
    est <- as.vector( my + xnew %*% be )  ## predicted values for PCA model
  }
  rownames(be) <- colnames(x)
  list(be = be, per = per[k], est = est)
}
