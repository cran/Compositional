################################
#### Multivariate or univariate regression with compositional data
#### in the covariates side using the alpha-transformation
#### Tuning the values of a and the number of principal components
#### Tsagris Michail 8/2015
#### mtsagris@yahoo.gr
#### References:Tsagris M. (2015)
#### Regression analysis with compositional data containing zero values
#### Chilean journal of statistics 6(2): 47-57
################################

alfapcr.tune <- function(y, x, M = 10, maxk = 50, a = seq(-1, 1, by = 0.1),
  oiko = "normal", seed = FALSE, ncores = 2) {
  ## oiko can be either "normal", "binomial" or "poisson"
  ## depending on the type of the independent variable
  ## "normal" is set by default
  x <- as.matrix(x)
  x <- x/ rowSums(x)
  d <- ncol(x) - 1
  if ( min(x) == 0 )  a <- a[a>0]  ## checks for zero values in the data.
  mspe <- matrix(nrow = length(a), ncol = d)
  if (oiko == 'normal') {
    for ( i in 1:length(a) ) {
      z <- alfa(x, a[i])$aff
      mod <- pcr.tune(y, z, M = M, maxk = maxk, seed = seed,
      ncores = ncores)
      mspe[i, ] <- mod$mspe
    }
  } else {
    for ( i in 1:length(a) ) {
      z <- alfa(x, a[i])$aff
      mod <- glmpcr.tune(y, z, M = M, maxk = maxk, oiko = oiko,
      seed = TRUE, ncores = ncores)
      mspe[i, ] <- mod$mspe
    }
  }
  rownames(mspe) <- a
  colnames(mspe) <- paste("PC", 1:d, sep = "")
  mspe
}
