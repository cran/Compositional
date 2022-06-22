################################
#### Principal components regression for binary and poisson regression
#### Selection of the number of principal components
#### via K-fold cross validation
#### Tsagris Michail 1/2016
#### mtsagris@yahoo.gr
################################
glmpcr.tune <- function(y, x, nfolds = 10, maxk = 10, folds = NULL, ncores = 1, seed = NULL, graph = TRUE) {
  ## y is the UNIVARIATE dependent variable
  ## y is either a binary variable (binary logistic regression)
  ## or a discrete variable (Poisson regression)
  ## x contains the independent variables
  ## fraction denotes the percentage of observations
  ## to be used as the test set
  ## the 1-fraction proportion of the data will be the training set
  ## R is the number of cross validations
  ## if ncores==1, then 1 processor is used, otherwise more are
  ## used (parallel computing)
  n <- dim(x)[1]
  p <- dim(x)[2]
  if ( maxk > p ) maxk <- p  ## just a check
  if ( is.null(folds) )  folds <- Compositional::makefolds(y, nfolds = nfolds,
                                                           stratified = FALSE, seed = seed)
  nfolds <- length(folds)
  msp <- matrix( nrow = nfolds, ncol = maxk )
  ## deigma will contain the positions of the test set
  ## this is stored but not showed in the end
  ## the user can access it though by running
  ## the commands outside this function
  if ( length( Rfast::sort_unique(y) ) == 2 ) {
    oiko <- "binomial"
  } else oiko <- "poisson"

  if (ncores <= 1) {
    runtime <- proc.time()
    for (vim in 1:nfolds) {
      ytest <- y[ folds[[ vim ]] ]   ## test set dependent vars
      ytrain <- y[ -folds[[ vim ]] ]   ## train set dependent vars
      xtrain <- x[ -folds[[ vim ]], , drop = FALSE]   ## train set independent vars
      xtest <- x[ folds[[ vim ]], , drop = FALSE ]  ## test set independent vars
	    vec <- prcomp(xtrain, center = FALSE)$rotation
      z <- xtrain %*% vec  ## PCA scores

      for ( j in 1:maxk) {
        if (oiko == "binomial") {
          be <- Rfast::glm_logistic(z[, 1:j], ytrain)$be
        } else {
          be <- Rfast::glm_poisson(z[, 1:j], ytrain)$be
        }
        ztest <- xtest %*% vec[, 1:j, drop = FALSE]  ## PCA scores
        es <- as.vector( ztest %*% be[-1] ) + be[1]

        if (oiko == "binomial") {
          est <- as.vector(  exp(es) / ( 1 + exp(es) )  )
          ri <-  -2 *( ytest * log(est) + (1 - ytest) * log(1 - est) )
        } else {
          est <- as.vector( exp(es) )
          ri <- 2 * ytest * log(ytest / est)
        }
        msp[vim, j] <- sum( ri, na.rm = TRUE )
      }
    }
    runtime <- proc.time() - runtime

  } else {
    runtime <- proc.time()
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    er <- numeric(maxk)
    if ( is.null(folds) )  folds <- Compositional::makefolds(y, nfolds = nfolds,
                                                             stratified = FALSE, seed = seed)
    msp <- foreach::foreach(vim = 1:nfolds, .combine = rbind, .packages = "Rfast", .export = c("glm_logistic", "glm_poisson") ) %dopar% {
      ytest <- y[ folds[[ vim ]] ]  ## test set dependent vars
      ytrain <-  y[ -folds[[ vim ]] ]   ## train set dependent vars
      xtrain <- x[ -folds[[ vim ]], , drop = FALSE]   ## train set independent vars
      xtest <- x[ folds[[ vim ]], , drop = FALSE]  ## test set independent vars
	  vec <- prcomp(xtrain, center = FALSE)$rotation
      z <- xtrain %*% vec  ## PCA scores

      for ( j in 1:maxk) {
        if (oiko == "binomial") {
          be <- Rfast::glm_logistic(z[, 1:j], ytrain)$be
        } else {
          be <- Rfast::glm_poisson(z[, 1:j], ytrain)$be
        }
        ztest <- xtest %*% vec[, 1:j, drop = FALSE]  ## PCA scores
        es <- as.vector( ztest %*% be[-1] ) + be[1]

        if (oiko == "binomial") {
          est <- exp(es) / ( 1 + exp(es) )
          ri <-  -2 *( ytest * log(est) + (1 - ytest) * log(1 - est) )
        } else {
          est <- exp(es)
          ri <- 2 * ytest * log(ytest / est)
        }
        er[j] <- sum( ri, na.rm = TRUE )
      }
      return(er)
    }
    stopCluster(cl)
    runtime <- proc.time() - runtime
  }

  mpd <- Rfast::colmeans(msp)
  if ( graph )  plot( 1:maxk, mpd, xlab = "Number of principal components",
                      ylab = "Mean predicted deviance", type = "b", cex.lab = 1.2,
                      cex.axis = 1.2, col = "green", pch = 16 )
  abline(v = 1:maxk, col = "lightgrey", lty = 2)
  abline(h = seq(min(mpd), max(mpd), length = 10), col = "lightgrey", lty = 2)

  names(mpd) <- paste("PC", 1:maxk, sep = " ")
  performance <- min(mpd)
  names(performance) <- "MPD"
  list(msp = msp, mpd = mpd, k = which.min(mpd), performance = performance, runtime = runtime)
}
