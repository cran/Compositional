################################
#### Principal components regression for binary and poisson regression
#### Selection of the number of principal components
#### via K-fold cross validation
#### Tsagris Michail 1/2016
#### mtsagris@yahoo.gr
################################
multinompcr.tune <- function(y, x, M = 10, maxk = 10, mat = NULL, ncores = 1, graph = TRUE) {
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
  if ( !is.numeric(y) )   y <- as.numeric(y)
  n <- dim(x)[1]
  p <- dim(x)[2]
  if ( maxk > p ) maxk <- p  ## just a check

  if ( is.null(mat) ) {
    nu <- sample(1:n, min( n, round(n / M) * M ) )
    ## It may be the case this new nu is not exactly the same
    ## as the one specified by the user
    ## to a matrix a warning message should appear
    options(warn = -1)
    mat <- matrix( nu, ncol = M )
  } else  mat <- mat

  M <- dim(mat)[2]
  rmat <- dim(mat)[1]
  msp <- matrix( nrow = M, ncol = maxk )

  if (ncores == 1) {
    runtime <- proc.time()
    for (vim in 1:M) {
      ytest <- y[mat[, vim] ]   ## test set dependent vars
      ytrain <- y[-mat[, vim] ]   ## train set dependent vars
      xtrain <- x[-mat[, vim], , drop = FALSE]   ## train set independent vars
      xtest <- x[mat[, vim], , drop = FALSE ]  ## test set independent vars
      vec <- prcomp(xtrain, center = FALSE)$rotation
      z <- xtrain %*% vec  ## PCA scores

      for ( j in 1:maxk) {
        mod <-try( Rfast::multinom.reg(ytrain, z[, 1:j]), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          est <- NULL
        } else {
          be <- mod$be
          ztest <- cbind(1, xtest %*% vec[, 1:j, drop = FALSE])  ## PCA scores
          es <- cbind(1, exp( ztest %*% be ) )
          est <- es / Rfast::rowsums(es)
          est <- Rfast::rowMaxs(est)
        }
        msp[vim, j] <- sum( est == ytest ) / rmat
      }
    }
    runtime <- proc.time() - runtime

  } else {
    runtime <- proc.time()
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    er <- numeric(maxk)
    msp <- foreach(vim = 1:M, .combine = rbind, .packages = "Rfast", .export = c("multinom.reg", "rowMaxs") ) %dopar% {
      ytest <- y[mat[, vim] ]  ## test set dependent vars
      ytrain <-  y[-mat[, vim] ]   ## train set dependent vars
      xtrain <- x[-mat[, vim], , drop = FALSE]   ## train set independent vars
      xtest <- x[mat[, vim], , drop = FALSE]  ## test set independent vars
      vec <- prcomp(xtrain, center = FALSE)$rotation
      z <- xtrain %*% vec  ## PCA scores

      for ( j in 1:maxk) {
        mod <-try( Rfast::multinom.reg(ytrain, z[, 1:j]), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          est <- NULL
        } else {
          be <- mod$be
          ztest <- cbind(1, xtest %*% vec[, 1:j, drop = FALSE])  ## PCA scores
          es <- cbind(1, exp( ztest %*% be ) )
          est <- es / Rfast::rowsums(es)
          est <- Rfast::rowMaxs(est)
        }
        er[j] <- sum( est == ytest ) / rmat
      }
      return(er)
    }
    stopCluster(cl)
    runtime <- proc.time() - runtime
  }

  mpd <- Rfast::colmeans(msp)
  if ( graph )  plot(1:maxk, mpd, xlab = "Number of principal components", ylab = "Mean predicted deviance", type = "b" , cex.lab = 1.3)

  names(mpd) <- paste("PC", 1:maxk, sep = " ")
  performance <- max(mpd)
  names(performance) <- "MPD"
  list(msp = msp, mpd = mpd, k = which.max(mpd), performance = performance, runtime = runtime)
}
