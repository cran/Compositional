################################
#### Multivariate or univariate regression with compositional data
#### in the covariates side using the alpha-transformation
#### Tuning the values of a and the number of principal components
#### Tsagris Michail 2/2016
#### mtsagris@yahoo.gr
#### References:Tsagris M. (2015)
#### Regression analysis with compositional data containing zero values
#### Chilean journal of statistics 6(2): 47-57
################################
alfaridge.tune <- function(y, x, nfolds = 10, a = seq(-1, 1, by = 0.1), lambda = seq(0, 2, by = 0.1), folds = NULL,
                           ncores = 1, graph = TRUE, col.nu = 15, seed = FALSE) {

  if ( min(x) == 0 )  a <- a[a>0]  ## checks for zero values in the data.
  da <- length(a)
  n <- dim(x)[1]
  ina <- 1:n
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                                           stratified = FALSE, seed = seed)
  nfolds <- length(folds)
  mspe2 <- array( dim = c( nfolds, length(lambda), da ) )

  tac <- proc.time()
  for ( i in 1:da ) {
    z <- alfa(x, a[i])$aff
    mod <- Compositional::ridge.tune(y, z, nfolds = nfolds, lambda = lambda, folds = folds, ncores = ncores,
                                     seed = seed, graph = FALSE)
    mspe2[, , i] <- mod$msp
  }

  dimnames(mspe2) <- list(folds = 1:nfolds, lambda = lambda, a = a)
  mspe <- array( dim = c(da, length(lambda), nfolds) )
  for (i in 1:nfolds)  mspe[, , i] <- t( mspe2[i, , 1:da] )
  dimnames(mspe) <- list(a = a, lambda = lambda, folds = 1:nfolds )
  mean.mspe <- apply(mspe, 1:2, mean)
  best.par <- ( which(mean.mspe == min(mean.mspe), arr.ind = TRUE)[1, ] )
  opt.mspe <- mean.mspe[ best.par[1], best.par[2] ]
  rownames(mean.mspe) <- a
  colnames(mspe) <- lambda
  if ( graph )  filled.contour( a, lambda, mean.mspe, xlab = expression( paste(alpha, " values") ),
                                ylab = expression( paste(lambda, " values") ), cex.lab = 1.3 )

  best.par <- c( a[ best.par[1] ], best.par[2] )
  names(best.par) <- c("alpha", "lambda")
  performance <- opt.mspe
  names(performance) <- "mspe"
  runtime <- proc.time()- tac
  list(mspe = mean.mspe, best.par = best.par, performance = performance, runtime = runtime)
}
