################################
#### Ridge regression tuning of the
#### lambda parameter via k-fold cross validation
#### Tsagris Michail 1/2016
#### mtsagris@yahoo.gr
################################

ridge.tune <- function(y, x, M = 10, lambda = seq(0, 2, by = 0.1),
                       seed = FALSE, ncores = 2, graph = TRUE) {
  ## y is the univariate dependent real valued variable
  ## x contains the independent variables(s)
  ## M is the number of folds, set to 10 by default
  ## lambda is a vector with a grid of values of lambda
  ## ncores is the number of cores to use
  y <- as.vector(y)  ## makes sure y is a matrix
  x <- as.matrix(x)
  n <- length(y)  ## sample size
  k <- length(lambda)
  ## if seed==TRUE then the results will always be the same
  if (seed == TRUE)  set.seed(1234567)
  nu <- sample(1:n, min( n, round(n / M) * M ) )
  ## It may be the case this new nu is not exactly the same
  ## as the one specified by the user
  options(warn = -1)
  mat <- matrix( nu, ncol = M ) # if the length of nu does not fit
  ## to a matrix a warning message should appear
  msp <- matrix( nrow = M, ncol = k )
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
      for ( j in 1:k ) {
        est <- ridge.reg(ytrain, xtrain, lambda[j], B = 1, xtest)$est
        msp[vim, j] <- sum( (ytest - est)^2 ) / rmat
      }
    }
  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    er <- numeric(k)
    msp <- foreach::foreach(vim = 1:M, .combine = rbind, .export = "ridge.reg") %dopar% {
      if (seed == TRUE)  set.seed( vim:c(vim + 10) ) ## If seed is TRUE the results
      ## will always be the same
      ytest <- as.vector( y[mat[, vim] ] )  ## test set dependent vars
      ytrain <- as.vector( y[-mat[, vim] ] )  ## train set dependent vars
      xtrain <- as.matrix( x[-mat[, vim], ] )  ## train set independent vars
      xtest <- as.matrix( x[mat[, vim], ] )  ## test set independent vars
      for ( j in 1:k ) {
        est <- ridge.reg(ytrain, xtrain, lambda[j], B = 1, xtest)$est
        er[j] <- sum( (ytest - est)^2 ) / rmat
      }
      return(er)
    }
    stopCluster(cl)
  }
  mspe <- colMeans(msp)
  bias <- msp[, which.min(mspe)] - apply(msp, 1, min)  ## TT estimate of bias
  estb <- mean( bias )  ## TT estimate of bias
  if (graph == TRUE) {
    plot(lambda, mspe, type = 'b', ylim = c( min(mspe), max(mspe) ),
    ylab = "Mean squared error of prediction",
    xlab = expression(paste(lambda, " values")) )
  }
  names(mspe) <- lambda
  rownames(msp) <- 1:M
  colnames(msp) <- lambda
  performance <- c( min(mspe) + estb, estb)
  names(performance) <- c("MSPE", "Estimated bias")
  list(msp = msp, mspe = mspe, lambda = which.min(mspe), performance = performance)
}
