################################
#### Principal components regression  for binary and poisson regression
#### Tsagris Michail 8/2015
#### mtsagris@yahoo.gr
################################
glm.pcr <- function(y, x, k = 1, xnew = NULL) {
  ## y is either a binary variable 0, 1 (binomial) or
  ## a numerical variable with counts (poisson)
  ## x contains the independent variables
  ## k shows the number of components to keep
  ## oiko can be "binomial" or "poisson"
  p <- dim(x)[2]
  m <- Rfast::colmeans(x)
  x <- Rfast::standardise(x)  ## standardize the independent variables
  ## eig <- eigen(crossprod(x))  ## eigen analysis of the design matrix
  eig <- prcomp(x, center = FALSE)
  values <- eig$sdev^2  ## eigenvalues
  per <- cumsum( values / sum(values) )  ## cumulative proportion of eigenvalues
  vec <- eig$rotation  ## eigenvectors, or principal components
  z <- x %*% vec  ## PCA scores

  if ( length( Rfast::sort_unique(y) ) == 2 ) {
    oiko <- "binomial"
    mod <- Rfast::glm_logistic(z[, 1:k], y, full = TRUE)
    b <- mod$info[, 1]
  } else {
    oiko <- "poisson"
    mod <- Rfast::glm_poisson(z[, 1:k], y, full = TRUE)
    b <- mod$info[, 1]
  }
  be <- vec[, 1:k, drop = FALSE] %*% b[-1]

  if ( !is.null(xnew) ) {
    xnew <- matrix(xnew, ncol = p)
    s <- Rfast::colVars(x, std = TRUE)
    xnew <- t( ( t(xnew) - m ) / s )## standardize the xnew values
    es <- as.vector( xnew %*% be ) + b[1]
  } else  es <- as.vector( x %*% be ) + b[1]

  if (oiko == "binomial") {
    est <- exp(es) / (1 + exp(es))
  } else est <- exp(es)     ## fitted values for PCA model

  list(model = mod, per = per[k], est = est)
}
