################################
#### Principal components regression  for binary and poisson regression
#### Tsagris Michail 8/2015
#### mtsagris@yahoo.gr
################################

glm.pcr <- function(y, x, k, oiko = "binomial", xnew = NULL) {
  ## y is either a binary variable 0, 1 (binomial) or 
  ## a numerical variable with counts (poisson)
  ## x contains the independent variables
  ## k shows the number of components to keep
  ## oiko can be "binomial" or "poisson"
  x <- as.matrix(x)
  y <- as.vector(y)
  n <- nrow(x)
  p <- ncol(x)
  m <- colMeans(x)
  s <- apply(x, 2, sd)
  s <- diag(1/s)
  x <- scale(x)[1:n, ]  ## standardize the independent variables
  eig <- eigen(crossprod(x))  ## eigen analysis of the design matrix
  values <- eig$values  ## eigenvalues
  per <- cumsum( values / sum(values) )  ## cumulative proportion of eigenvalues
  vec <- eig$vectors  ## eigenvectors, or principal components
  z <- x %*% vec  ## PCA scores
  mod <- glm(y ~ z[ , 1:k], family = oiko )
  if ( !is.null(xnew) ) {
    xnew <- as.matrix(xnew)
    xnew <- matrix(xnew, ncol = p)
    nu <- nrow(xnew)
    xnew <- ( xnew - rep(m, rep(nu, p)) ) %*% s ## standardize the xnew values 
    znew <- cbind(1, xnew %*% vec) 
    znew <- znew[, 1:c(k + 1)]
    be <- coef(mod)
    es <- znew %*% be 
    if (oiko == "binomial") {
      est <- as.vector(  exp(es)/(1 + exp(es))  )
    } else est <- as.vector( exp(es) )
  } else {
    est <- as.vector( fitted(mod) )  ## fitted values for PCA model
  }
  list(model = summary(mod), per = per[k], est = est)
}