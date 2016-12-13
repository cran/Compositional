################################
#### Ridge regression tuning of the
#### lambda parameter via k-fold cross validation
#### Tsagris Michail 1/2016
#### mtsagris@yahoo.gr
################################

ridge.tune <- function(y, x, M = 10, lambda = seq(0, 2, by = 0.1),
                       mat = NULL, ncores = 1, graph = FALSE) {

  ## y is the univariate or multivariate dependent variable
  ## x contains the independent variables(s)
  ## M is the number of folds, set to 10 by default
  ## lambda is a vector with a grid of values of lambda
  ## ncores is the number of cores to use

  y <- as.matrix(y)
  n <- dim(y)[1]  ## sample size
  k <- length(lambda)
  di <- dim(y)[2]  ## dimensionality of y
  p <- dim(x)[2]  ## dimensionality of x
  if ( is.null(mat) ) {
    nu <- sample(1:n, min( n, round(n / M) * M ) )
    ## It may be the case this new nu is not exactly the same
    ## as the one specified by the user
    ## to a matrix a warning message should appear
    options(warn = -1)
    mat <- matrix( nu, ncol = M ) # if the length of nu does not fit
  } else  mat <- mat

  M <- dim(mat)[2]
  rmat <- dim(mat)[1]
  msp <- matrix( nrow = M, ncol = k)
  ## deigma will contain the positions of the test set
  ## this is stored but not showed in the end
  ## the user can access it though by running
  ## the commands outside this function
  if (ncores == 1) {
    runtime <- proc.time()
    for (vim in 1:M) {
      ytest <- y[ mat[, vim] ]   ## test set dependent vars
      ytrain <- y[ -mat[, vim], ]   ## train set dependent vars
      my <- sum(ytrain) / rmat
      yy <- ytrain - my  ## center the dependent variables

      xtrain <- as.matrix( x[ -mat[, vim], ] )  ## train set independent vars
      mx <- Rfast::colmeans(xtrain)
      xtest <- as.matrix( x[ mat[, vim], ] )  ## test set independent vars
      s <- Rfast::colVars(xtrain, std = TRUE)
      xtest <- t( ( t(xtest) - mx ) / s ) ## standardize the xtest
      xx <- t( ( t(xtrain) - mx ) / s ) ## standardize the independent variables

      sa <- svd(xx)
      d <- sa$d    ;    v <- t(sa$v)    ;     tu <- t(sa$u)
      d2 <- d^2    ;    A <- d * tu %*% yy
      for (i in 1:k) {
        ## beta <- ( v %*% (tu * ( d / ( d^2 + lambda[i] ) ) ) ) %*% yy
        beta <- crossprod( v / ( d2 + lambda[i] ), A )
        est <- xtest %*% beta + my
        msp[vim, i] <- sum( (ytest - est)^2 ) / rmat
      }
    }
    runtime <- proc.time() - runtime

  } else {
    runtime <- proc.time()
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    pe <- numeric(k)

    msp <- foreach(vim = 1:M, .combine = rbind, .packages = "Rfast", .export = c("colmeans", "colVars") ) %dopar% {
      ytest <- y[ mat[, vim], ]   ## test set dependent vars
      ytrain <- y[ -mat[, vim], ]   ## train set dependent vars
      my <- sum(ytrain) / rmat
      yy <- ytrain- my  ## center the dependent variables

      xtrain <- as.matrix( x[ -mat[, vim], ] )  ## train set independent vars
      mx <- Rfast::colmeans(xtrain)
      xtest <- as.matrix( x[ mat[, vim], ] )  ## test set independent vars
      s <- Rfast::colVars(xtrain, std = TRUE)
      xtest <- t( ( t(xtest) - mx ) / s ) ## standardize the xtest
      xx <- t( ( t(xtrain) - mx ) / s ) ## standardize the independent variables

      sa <- svd(xx)
      d <- sa$d    ;    v <- t(sa$v)    ;     tu <- t(sa$u)
      d2 <- d^2    ;    A <- d * tu %*% yy
      for ( i in 1:k ) {
        ## beta <- ( v %*% (tu * ( d / ( d^2 + lambda[i] ) ) ) ) %*% yy
        beta <- crossprod( v / ( d2 + lambda[i] ), A )
        est <- xtest %*% beta + my
        pe[i] <- sum( (ytest - est)^2 ) / rmat
      }
      return(pe)
    }
    stopCluster(cl)
    runtime <- proc.time() - runtime
  }

  mspe <- Rfast::colmeans(msp)
  bias <- msp[ , which.min(mspe)] - Rfast::rowMins(msp, value = TRUE)   ## apply(msp, 1, min)  ## TT estimate of bias
  estb <- mean( bias )  ## TT estimate of bias

  if ( graph ) {
    plot(lambda, mspe, type = 'b', ylim = c(min(mspe), max(mspe)),
         ylab = "Mean squared error of prediction",
         xlab = expression(paste(lambda, " values")) )
  }

  names(mspe) <- lambda
  performance <- c( min(mspe) + estb, estb)
  names(performance) <- c("MSPE", "Estimated bias")
  list(msp = msp, mspe = mspe, lambda = which.min(mspe), performance = performance,
       runtime = runtime)

}
