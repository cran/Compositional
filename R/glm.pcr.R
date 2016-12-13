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

  n <- length(y)
  p <- dim(x)[2]
  m <- Rfast::colmeans(x)

  x <- Rfast::standardise(x)  ## standardize the independent variables
  eig <- eigen(crossprod(x))  ## eigen analysis of the design matrix
  values <- eig$values  ## eigenvalues
  per <- cumsum( values / sum(values) )  ## cumulative proportion of eigenvalues
  vec <- eig$vectors  ## eigenvectors, or principal components
  z <- x %*% vec  ## PCA scores

  if ( length( Rfast::sort_unique(y) ) == 2 ) {
    oiko <- "binomial"
  } else oiko <- "poisson"

  mod <- glm(y ~ z[, 1:k], family = oiko )
  b <- coef(mod)
  be <- vec[, 1:k] %*% as.matrix( b[-1] )

  if ( !is.null(xnew) ) {
    xnew <- matrix(xnew, ncol = p)
    s <- Rfast::colVars(x, std = TRUE)
    xnew <- t( ( t(xnew) - m ) / s )## standardize the xnew values
    es <- as.vector( xnew %*% be ) + b[1]
  } else  es <- as.vector( x %*% be ) + b[1]

  if (oiko == "binomial") {
    est <- exp(es) / (1 + exp(es)) 
  } else est <- exp(es)     ## fitted values for PCA model

  list(model = summary(mod), per = per[k], est = est)
}
