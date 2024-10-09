################################
#### Multivariate or univariate regression with compositional data
#### in the covariates side using the alpha-transformation
#### Tsagris Michail 8/2015
#### mtsagris@yahoo.gr
#### References:Tsagris M. (2015)
#### Regression analysis with compositional data containing zero values
#### Chilean journal of statistics 6(2): 47-57
################################
alfa.pcr <- function(y, x, a, k, model = "gaussian", xnew = NULL) {
  ## y is dependent univariate variable. It can be matrix or vector
  ## x are compositional data, the covariates
  ## a is the value of a for the alpha-transformation
  ## k is the number of principal components to use
  ## oiko can be either "normal", "binomial" or "poisson"
  ## depending on the type of the independent variable
  ## "normal" is set by default
  x <- Compositional::alfa(x, a, h = TRUE)$aff ## apply the alpha-transformation
  dm <- dim(x)
  p <- dm[2]
  if ( length(k) == 1  & k > p ) k <- p
  if ( length(k)> 1  &  max(k) > p )  k <- 1:p

  if ( model == 'gaussian' ) {
    if ( !is.null(xnew) ) {
      xnew <- Compositional::alfa(xnew, a, h = TRUE)$aff ## apply the alpha-transformation
    }
    mod <- Rfast2::pcr(y, x, k, xnew = xnew)

  } else if ( model == "multinomial" ) {
    p <- dim(x)[2]
    eig <- prcomp(x, center = FALSE, scale = FALSE)
    values <- eig$sdev^2  ## eigenvalues
    per <- cumsum( values / sum(values) )  ## cumulative proportion of eigenvalues
    vec <- eig$rotation[, 1:k, drop = FALSE]  ## eigenvectors, or principal components
    z <- eig$x[, 1:k, drop = FALSE]  ## PCA scores
    mod <- try( Rfast::multinom.reg(y, z), silent = TRUE)
    if ( identical( class(mod), "try-error" ) ) {
      be <- NULL
      est <- NULL
    } else {
      be <- mod$be
      est <- NULL
      if ( !is.null(xnew) ) {
        xnew <- matrix(xnew, ncol = p + 1)
        xnew <- Compositional::alfa(xnew, a, h = TRUE)$aff ## apply the alpha-transformation
        znew <- cbind(1, xnew %*% vec)  ## PCA scores
        es <- cbind(1, exp( znew %*% be ) )
        est <- es / Rfast::rowsums(es)
      }
    }
    mod <- list(model = mod, per = per[k], vec = vec, est = est)

  } else if ( model == "poisson" | model == "binomial" ) {
    if ( !is.null(xnew) ) {
      xnew <- matrix(xnew, ncol = p + 1)
      xnew <- Compositional::alfa(xnew, a, h = TRUE)$aff ## apply the alpha-transformation
    }
    mod <- Compositional::glm.pcr(y, x, k, xnew = xnew)
  }

  mod
}
