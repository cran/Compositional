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
  eig <- prcomp(x, center = FALSE)
  values <- eig$sdev^2  ## eigenvalues
  per <- cumsum( values / sum(values) )  ## cumulative proportion of eigenvalues
  vec <- eig$rotation[, 1:k, drop = FALSE]  ## eigenvectors, or principal components
  if ( !is.matrix(vec) )  vec <- as.matrix(vec)
  z <- x %*% vec  ## PCA scores

  if ( length( Rfast::sort_unique(y) ) == 2 ) {
    oiko <- "binomial"
    mod <- Rfast::glm_logistic(z, y)
    be <- mod$be
  } else {
    oiko <- "poisson"
    mod <- Rfast::glm_poisson(z, y)
    be <- mod$be
  }

  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- matrix(xnew, ncol = p)
    znew <- xnew %*% vec  ## PCA scores
    es <- as.vector( znew %*% be[-1] ) + be[1]
	if (oiko == "binomial") {
      est <- exp(es) / (1 + exp(es))
    } else est <- exp(es)     ## fitted values for PCA model
  }
  nam <- colnames(x)
  if ( is.null(nam) )  nam <- paste("X", 1:p, sep = "")
  rownames(be) <- c("constant", paste("PC", 1:k, sep = "") )
  rownames(vec) <- nam
  colnames(vec) <- paste("PC", 1:k, sep = "")
  list(model = mod, per = per[k], vec = vec, est = est)
}
