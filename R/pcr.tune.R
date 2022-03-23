################################
#### Selection of the number of principal components in PCR
#### via K-fold cross validation
#### Tsagris Michail 12/2013
#### mtsagris@yahoo.gr
#### References: Jolliffe I.T. (2002)
#### Principal Component Analysis p. 167-188.
################################
pcr.tune <- function(y, x, nfolds = 10, maxk = 50, folds = NULL, ncores = 1, seed = NULL, graph = TRUE) {
  ## y is the univariate dependent variable
  ## x contains the independent variables(s)
  ## M is the number of folds, set to 10 by default
  ## maxk is the maximum number of eigenvectors to conside
  ## ncores specifies how many cores to use
  n <- length(y)  ## sample size
  p <- dim(x)[2]  ## number of independent variables
  if ( maxk > p )  maxk <- p  ## just a check
  if ( is.null(folds) )  folds <- Compositional::makefolds(y, nfolds = nfolds,
                                                           stratified = FALSE, seed = seed)
  nfolds <- length(folds)

  if (ncores <= 1) {

    runtime <- proc.time()
    msp <- matrix( nrow = nfolds, ncol = maxk )

    for (vim in 1:nfolds) {
      ytest <- y[ folds[[ vim ]] ]  ## test set dependent vars
      ytrain <- y[ -folds[[ vim ]] ]   ## train set dependent vars
      xtrain <- x[ -folds[[ vim ]], , drop = FALSE]   ## train set independent vars
      xtest <- x[ folds[[ vim ]], , drop = FALSE]  ## test set independent vars
      est <- Compositional::pcr(ytrain, xtrain, k = 1:maxk, xnew = xtest)$est
      msp[vim, ] <- Rfast::colmeans( (est - ytest)^2 )
    }

    runtime <- proc.time() - runtime

  } else {

    runtime <- proc.time()

    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    er <- numeric(maxk)
    if ( is.null(folds) )  folds <- Compositional::makefolds(y, nfolds = nfolds,
                                                             stratified = FALSE, seed = seed)
    msp <- foreach::foreach(vim = 1:nfolds, .combine = rbind, .packages = c("Rfast", "Compositional") ) %dopar% {
      ytest <-  y[ folds[[ vim ]] ]  ## test set dependent vars
      ytrain <- y[ -folds[[ vim ]] ]   ## train set dependent vars
      xtrain <- x[ -folds[[ vim ]], , drop = FALSE]   ## train set independent vars
      xtest <- x[ folds[[ vim ]], , drop = FALSE]  ## test set independent vars
      est <- Compositional::pcr(ytrain, xtrain, k = 1:maxk, xnew = xtest)$est
      er <- Rfast::colmeans( (est - ytest)^2 )
      return(er)
    }
    parallel::stopCluster(cl)

    runtime <- proc.time() - runtime
  }

  mspe <- Rfast::colmeans(msp)
  if ( graph )  plot(1:maxk, mspe, xlab = "Number of principal components", ylab = "MSPE", type = "b",
                     cex.lab = 1.2, cex.axis = 1.2, col = "green", pch = 16)
  abline(v = 1:maxk, col = "lightgrey", lty = 2)
  abline(h = seq(min(mspe), max(mspe), length = 10), col = "lightgrey", lty = 2)

  names(mspe) <- paste("PC", 1:maxk, sep = " ")
  performance <- min(mspe)
  names(performance) <- "MSPE"
  list(msp = msp, mspe = mspe, k = which.min(mspe), performance = performance, runtime = runtime)
}
