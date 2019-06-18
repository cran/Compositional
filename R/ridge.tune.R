################################
#### Ridge regression tuning of the
#### lambda parameter via k-fold cross validation
#### Tsagris Michail 1/2016
#### mtsagris@yahoo.gr
################################
ridge.tune <- function(y, x, nfolds = 10, lambda = seq(0, 2, by = 0.1), folds = NULL, ncores = 1, seed = FALSE, graph = FALSE) {
  ## x contains the independent variables(s)
  ## nfolds is the number of folds, set to 10 by default
  ## lambda is a vector with a grid of values of lambda
  ## ncores is the number of cores to use
  n <- length(y)  ## sample size
  k <- length(lambda)
  di <- dim(y)[2]  ## dimensionality of y
  p <- dim(x)[2]  ## dimensionality of x
  ina <- 1:n
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                                           stratified = FALSE, seed = seed)
  nfolds <- length(folds)
  msp <- matrix( nrow = nfolds, ncol = k)
  ## deigma will contain the positions of the test set
  ## this is stored but not showed in the end
  ## the user can access it though by running
  ## the commands outside this function
  if (ncores <= 1) {
    runtime <- proc.time()
    for (vim in 1:nfolds) {
      ytest <- y[ folds[[ vim ]] ]   ## test set dependent vars
      ytrain <- y[ -folds[[ vim ]] ]   ## train set dependent vars
      my <- mean(ytrain)
      yy <- ytrain - my  ## center the dependent variables
      xtrain <- x[ -folds[[ vim ]], ]   ## train set independent vars
      xtest <- x[ folds[[ vim ]], , drop = FALSE]   ## test set independent vars
      sa <- svd(xtrain)
      d <- sa$d    ;    v <- t(sa$v)    ;     tu <- t(sa$u)
      d2 <- d^2    ;    A <- d * tu %*% yy
      for (i in 1:k) {
        ## beta <- ( v %*% (tu * ( d / ( d^2 + lambda[i] ) ) ) ) %*% yy
        beta <- crossprod( v / ( d2 + lambda[i] ), A )
        est <- xtest %*% beta + my
        msp[vim, i] <- mean( (ytest - est)^2 )
      }
    }
    runtime <- proc.time() - runtime

  } else {
    runtime <- proc.time()
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    pe <- numeric(k)
    if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                                             stratified = FALSE, seed = seed)
    nfolds <- length(folds)
    msp <- foreach::foreach(vim = 1:nfolds, .combine = rbind, .packages = "Rfast", .export = c("colmeans", "colVars") ) %dopar% {
      ytest <- y[ folds[[ vim ]] ]   ## test set dependent vars
      ytrain <- y[ -folds[[ vim ]] ]   ## train set dependent vars
      my <- mean(ytrain)
      yy <- ytrain - my  ## center the dependent variables
      xtrain <- x[ -folds[[ vim ]], ] ## train set independent vars
      xtest <- x[ folds[[ vim ]], , drop = FALSE]  ## test set independent vars
      sa <- svd(xtrain)
      d <- sa$d    ;    v <- t(sa$v)    ;     tu <- t(sa$u)
      d2 <- d^2    ;    A <- d * tu %*% yy
      for ( i in 1:k ) {
        ## beta <- ( v %*% (tu * ( d / ( d^2 + lambda[i] ) ) ) ) %*% yy
        beta <- crossprod( v / ( d2 + lambda[i] ), A )
        est <- xtest %*% beta + my
        pe[i] <- mean( (ytest - est)^2 )
      }
      return(pe)
    }
    stopCluster(cl)
    runtime <- proc.time() - runtime
  }

  mspe <- Rfast::colmeans(msp)
  if ( graph ) {
    plot( lambda, mspe, type = 'b', ylim = c( min(mspe), max(mspe) ),
         ylab = "Mean squared error of prediction", xlab = expression( paste(lambda, " values") ), cex.lab = 1.3 )
  }

  names(mspe) <- lambda
  performance <- min(mspe)
  names(performance) <- "MSPE"
  list(msp = msp, mspe = mspe, lambda = which.min(mspe), performance = performance, runtime = runtime)
}
