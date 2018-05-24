################################
#### Multivariate or univariate regression with compositional data
#### in the covariates side using the alpha-transformation
#### Tsagris Michail 8/2015
#### mtsagris@yahoo.gr
#### References:Tsagris M. (2015)
#### Regression analysis with compositional data containing zero values
#### Chilean journal of statistics 6(2): 47-57
################################
alfa.pcr <- function(y, x, a, k, xnew = NULL) {
  ## y is dependent univariate variable. It can be matrix or vector
  ## x are compositional data, the covariates
  ## a is the value of a for the alpha-transformation
  ## k is the number of principal components to use
  ## oiko can be either "normal", "binomial" or "poisson"
  ## depending on the type of the independent variable
  ## "normal" is set by default
  x <- alfa(x, a, h = TRUE)$aff ## apply the alpha-transformation
  dm <- dim(x)
  p <- dm[2]
  if (k > p)   k <- p

  if ( length( unique(y) ) == 2 ) {
    oiko <- "binomial"
  } else if ( sum(y - round(y) ) == 0 ) {
    oiko <- "poisson"
  } else oiko <- "normal"

  if (oiko == 'normal') {
    ## k shows the number of components to keep
    eig <- prcomp(x, center = FALSE, scale = FALSE)
    values <- eig$sdev^2
    per <- cumsum( values / sum(values) )  ## cumulative proportion of each eigenvalue
    vec <- eig$rotation[, 1:k, drop = FALSE]
    z <- cbind(1, x %*% vec)  ## PCA scores
    be <- solve( crossprod(z), crossprod(z, y) )
    est <- NULL
    if ( !is.null(xnew) ) {
      xnew <- matrix(xnew, ncol = p + 1)
      xnew <- alfa(xnew, a, h = TRUE)$aff ## apply the alpha-transformation
      znew <- xnew %*% vec  ## PCA scores
      est <- as.vector( znew %*% be[-1] ) + be[1]  ## predicted values for PCA model
    }
    rownames(be) <- colnames(x)
    mod <- list(be = be, per = per[k], vec = vec, est = est)


  } else {
    p <- dim(x)[2]
    eig <- prcomp(x, center = FALSE, scale = FALSE)
    values <- eig$sdev^2  ## eigenvalues
    per <- cumsum( values / sum(values) )  ## cumulative proportion of eigenvalues
    vec <- eig$rotation[, 1:k, drop = FALSE]  ## eigenvectors, or principal components
    z <- x %*% vec  ## PCA scores
    if ( length( Rfast::sort_unique(y) ) == 2 ) {
      oiko <- "binomial"
      mod <- Rfast::glm_logistic(z, y, full = TRUE)
      be <- mod$info[, 1]
    } else {
      oiko <- "poisson"
      mod <- Rfast::glm_poisson(z, y, full = TRUE)
      be <- mod$info[, 1]
    }
    est <- NULL
    if ( !is.null(xnew) ) {
      xnew <- matrix(xnew, ncol = p + 1)
      xnew <- alfa(xnew, a, h = TRUE)$aff ## apply the alpha-transformation
      znew <- xnew %*% vec  ## PCA scores
      es <- as.vector( znew %*% be[-1] ) + be[1]
      if (oiko == "binomial") {
        est <- exp(es) / (1 + exp(es))
      } else est <- exp(es)     ## fitted values for PCA model
    }
    mod <- list(model = mod, per = per[k], vec = vec, est = est)
  }

  mod
}
