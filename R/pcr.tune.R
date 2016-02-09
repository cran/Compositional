################################
#### Selection of the number of principal components in PCR
#### via K-fold cross validation
#### Tsagris Michail 12/2013
#### mtsagris@yahoo.gr
#### References: Jolliffe I.T. (2002)
#### Principal Component Analysis p. 167-188.
################################

pcr.tune <- function(y, x, M = 10, maxk = 50, seed = FALSE, ncores = 2) {
  ## y is the univariate dependent variable
  ## x contains the independent variables(s)
  ## M is the number of folds, set to 10 by default
  ## maxk is the maximum number of eigenvectors to conside
  ## ncores specifies how many cores to use
  y <- as.vector(y)  ## makes sure y is a vector
  x <- as.matrix(x)
  n <- length(y)  ## sample size
  p <- ncol(x)  ## number of independent variables
  if ( maxk > p ) maxk <- p  ## just a check
  ## if seed==TRUE then the results will always be the same
  if (seed == TRUE)  set.seed(1234567)
  nu <- sample(1:n, min( n, round(n / M) * M ) )
  ## It may be the case this new nu is not exactly the same
  ## as the one specified by the user
  options(warn = -1)
  mat <- matrix( nu, ncol = M ) # if the length of nu does not fit
  ## to a matrix a warning message should appear
  msp <- matrix( nrow = M, ncol = maxk )
  rmat <- nrow(mat)
  ## deigma will contain the positions of the test set
  ## this is stored but not showed in the end
  ## the user can access it though by running
  ## the commands outside this function
  if (ncores == 1) {
    for (vim in 1:M) {
      ytest <- as.vector( y[mat[, vim] ] )  ## test set dependent vars
      ytrain <- as.vector( y[-mat[, vim] ] )  ## train set dependent vars
      xtrain <- as.matrix( x[-mat[, vim], ] )  ## train set independent vars
      xtest <- as.matrix( x[mat[, vim], ] )  ## test set independent vars
      for ( j in 1:maxk ) {
        est <- pcr(ytrain, xtrain, j, xtest)$est
        msp[vim, j] <- sum( (ytest - est)^2 ) / rmat
      }
    }
  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    er <- numeric(maxk)
    msp <- foreach::foreach(vim = 1:M, .combine = rbind, .export = "pcr") %dopar% {
      if (seed == TRUE)  set.seed( vim:c(vim+10) ) ## If seed is TRUE the results
      ## will always be the same
      ytest <- as.vector( y[mat[, vim] ] )  ## test set dependent vars
      ytrain <- as.vector( y[-mat[, vim] ] )  ## train set dependent vars
      xtrain <- as.matrix( x[-mat[, vim], ] )  ## train set independent vars
      xtest <- as.matrix( x[mat[, vim], ] )  ## test set independent vars
      for ( j in 1:maxk ) {
        est <- pcr(ytrain, xtrain, j, xtest)$est
        er[j] <- sum( (ytest - est)^2 ) / rmat
      }
      return(er)
    }
    stopCluster(cl)
  }
  mspe <- colMeans(msp)
  bias <- msp[ ,which.min(mspe)] - apply(msp, 1, min)  ## TT estimate of bias
  estb <- mean( bias )  ## TT estimate of bias
  plot(1:maxk, mspe, xlab = "Number of principal components",
  ylab = "MSPE", type = "b")
  names(mspe) <- paste("PC", 1:maxk, sep = " ")
  performance <- c( min(mspe) + estb, estb)
  names(performance) <- c("MSPE", "Estimated bias")
  list(mspe = mspe, k = which.min(mspe), performance = performance)
}
