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
alfapcr.tune <- function(y, x, model = "gaussian", nfolds = 10, maxk = 50, a = seq(-1, 1, by = 0.1),
                         folds = NULL, ncores = 1, graph = TRUE, col.nu = 15, seed = FALSE) {
  ## model can be either "normal", "binomial" or "poisson"
  ## depending on the type of the independent variable
  ## "normal" is set by default
  n <- dim(x)[1]
  d <- dim(x)[2] - 1
  if ( min(x) == 0 )   a <- a[ a > 0 ]  ## checks for zero values in the data.
  da <- length(a)
  ina <- 1:n
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                                           stratified = FALSE, seed = seed)
  nfolds <- length(folds)
  mspe2 <- array( dim = c( nfolds, d, da) )

  if ( model == 'gaussian' ) {
    tic <- proc.time()
    for ( i in 1:da ) {
      z <- Compositional::alfa(x, a[i])$aff
      mod <- Compositional::pcr.tune(y, z, nfolds = nfolds, maxk = maxk, folds = folds, ncores = ncores, seed = seed, graph = FALSE)
      mspe2[, , i] <- mod$msp
    }
    toc <- proc.time() - tic

  } else if ( model == "multinomial" )  {
    tic <- proc.time()
    for ( i in 1:da ) {
      z <- Compositional::alfa(x, a[i])$aff
      mod <- Compositional::multinompcr.tune(y, z, nfolds = nfolds, maxk = maxk, folds = folds, ncores = ncores, seed = seed, graph = FALSE)
      mspe2[, , i] <- mod$msp
    }
    toc <- proc.time() - tic

  } else if ( model == "binomial"  |  model == "poisson" ) {
    tic <- proc.time()
    for ( i in 1:da ) {
      z <- Compositional::alfa(x, a[i])$aff
      mod <- Compositional::glmpcr.tune(y, z, nfolds = nfolds, maxk = maxk, folds = folds, ncores = ncores, seed = seed, graph = FALSE)
      mspe2[, , i] <- mod$msp
    }
    toc <- proc.time() - tic
  }

  dimnames(mspe2) <- list(folds = 1:nfolds, PC = paste("PC", 1:d, sep = ""), a = a)
  mspe <- array( dim = c(da, d, nfolds) )
  for (i in 1:nfolds)  mspe[, , i] <- t( mspe2[i, , 1:da] )
  dimnames(mspe) <- list(a = a, PC = paste("PC", 1:d, sep = ""), folds = 1:nfolds )
  mean.mspe <- t( colMeans( aperm(mspe) ) )   ## apply(mspe, 1:2, mean)
  if ( model == "multinomial" ) {
    best.par <- which(mean.mspe == max(mean.mspe), arr.ind = TRUE)[1, ]
  } else  best.par <- which(mean.mspe == min(mean.mspe), arr.ind = TRUE)[1, ]
  performance <- mean.mspe[ best.par[1], best.par[2] ]
  names(performance) <- "mspe"
  rownames(mean.mspe) <- a
  colnames(mspe) <- paste("PC", 1:d, sep = "")

  if ( graph )  filled.contour(a, 1:d, mean.mspe, xlab = expression( paste(alpha, " values") ),
                               ylab = "Number of PCs", cex.lab = 1.2, cex.axis = 1.2)
  best.par <- c( a[ best.par[1] ], best.par[2] )
  names(best.par) <- c("alpha", "PC")
  list(mspe = mean.mspe, best.par = best.par, performance = performance, runtime = toc)
}
