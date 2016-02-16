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
  oiko = "normal", seed = FALSE, ncores = 2, graph = TRUE, col.nu = 15) {
  ## oiko can be either "normal", "binomial" or "poisson"
  ## depending on the type of the independent variable
  ## "normal" is set by default
  x <- as.matrix(x)
  x <- x/ rowSums(x)
  d <- ncol(x) - 1
  if ( min(x) == 0 )  a <- a[a>0]  ## checks for zero values in the data.
  da <- length(a)
  mspe2 <- array( dim = c( M, d, da) )
  if (oiko == 'normal') {
    for ( i in 1:da ) {
      z <- alfa(x, a[i])$aff
      mod <- pcr.tune(y, z, M = M, maxk = maxk, seed = seed,
      ncores = ncores, graph = FALSE)
      mspe2[, , i] <- mod$msp
    }
  } else {
    for ( i in 1:da ) {
      z <- alfa(x, a[i])$aff
      mod <- glmpcr.tune(y, z, M = M, maxk = maxk, oiko = oiko,
      seed = TRUE, ncores = ncores, graph = FALSE)
      mspe2[, , i] <- mod$msp
    }
  }
  dimnames(mspe2) <- list(folds = 1:M, PC = paste("PC", 1:d, sep = ""), a = a)
  mspe <- array( dim = c(da, d, M) )
  for (i in 1:M)  mspe[, , i] <- t( mspe2[i, , 1:da] )
  dimnames(mspe) <- list(a = a, PC = paste("PC", 1:d, sep = ""), folds = 1:M )
  mean.mspe <- apply(mspe, 1:2, mean)
  best.par <- ( which(mean.mspe == min(mean.mspe), arr.ind = TRUE)[1, ] )
  opt.mspe <- mean.mspe[ best.par[1], best.par[2] ]
  estb <- mspe[ best.par[1], best.par[2], 1:M ] - apply(mspe, 3, min)
  bias <- mean(estb)
  rownames(mean.mspe) = a   ;  colnames(mspe) = paste("PC", 1:d, sep = "")
  if (graph == TRUE) {
    filled.contour(a, 1:d, mean.mspe, col = terrain.colors(col.nu),
     xlab = expression( paste(alpha, " values") ), ylab = "Number of PCs")
  }
  best.par <- c( a[ best.par[1] ], best.par[2] )
  names(best.par) <- c("alpha", "PC")
  performance <- c(opt.mspe, bias)
  names(performance) <- c("bias corrected mspe", "estimated bias")
  list(mspe = mean.mspe, best.par = best.par, performance = performance)
}
