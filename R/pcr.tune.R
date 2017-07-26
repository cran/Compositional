################################
#### Selection of the number of principal components in PCR
#### via K-fold cross validation
#### Tsagris Michail 12/2013
#### mtsagris@yahoo.gr
#### References: Jolliffe I.T. (2002)
#### Principal Component Analysis p. 167-188.
################################
pcr.tune <- function(y, x, M = 10, maxk = 50, mat = NULL, ncores = 1, graph = TRUE) {
  ## y is the univariate dependent variable
  ## x contains the independent variables(s)
  ## M is the number of folds, set to 10 by default
  ## maxk is the maximum number of eigenvectors to conside
  ## ncores specifies how many cores to use
  n <- length(y)  ## sample size
  p <- dim(x)[2]  ## number of independent variables
  if ( maxk > p )  maxk <- p  ## just a check

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
  ## deigma will contain the positions of the test set
  ## this is stored but not showed in the end
  ## the user can access it though by running
  ## the commands outside this function
  if (ncores == 1) {

    runtime <- proc.time()

    for (vim in 1:M) {

      ytest <- y[mat[, vim] ]  ## test set dependent vars
      ytrain <- y[-mat[, vim] ]   ## train set dependent vars
      xtrain <- x[-mat[, vim], ]   ## train set independent vars
      xtest <- x[mat[, vim], , drop = FALSE]  ## test set independent vars

      m <- mean(ytrain)
      ytrain <- ytrain - m  ## standardize the dependent variable
      mx <- Rfast::colmeans(xtrain)
      s <- Rfast::colVars(xtrain, std = TRUE)
      #mtrain <- t( xtrain )
      #mtrain <- mtrain - mx
      #mtrain <- t( mtrain / sqrt( Rfast::rowsums(mtrain^2) ) )
      #sar <- tcrossprod( mtrain )
      #eig <- eigen( sar )  ## eigen analysis of the design matrix
      #vec <- eig$vectors  ## eigenvectors, or principal components
	  vec <- prcomp(xtrain, center = TRUE, scale = TRUE)$rotation
      z <- xtrain %*% vec  ## PCA scores
      xnew <- t( ( t(xtest) - mx ) / s ) ## standardize the xnew values

      for ( j in 1:maxk ) {
        zzk <- crossprod(z[, 1:j])
        be <- vec[, 1:j] %*% solve( zzk, crossprod( z[, 1:j], ytrain ) )
        ## b is the PCA based coefficients
        est <- as.vector( m + xnew %*% be )  ## predicted values for PCA model
        msp[vim, j] <- sum( (ytest - est)^2 ) / rmat
      }
    }

    runtime <- proc.time() - runtime

  } else {

    runtime <- proc.time()

    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    er <- numeric(maxk)
    msp <- foreach::foreach(vim = 1:M, .combine = rbind, .packages = "Rfast",
                            .export = c("colVars", "colmeans") ) %dopar% {
      ytest <-  y[mat[, vim] ]  ## test set dependent vars
      ytrain <- y[-mat[, vim] ]   ## train set dependent vars
      xtrain <- x[-mat[, vim], ]   ## train set independent vars
      xtest <-  x[mat[, vim], , drop = FALSE]  ## test set independent vars
      m <- mean(ytrain)
      ytrain <- ytrain - m  ## standardize the dependent variable
      mx <- Rfast::colmeans(xtrain)
      s <- Rfast::colVars(xtrain, std = TRUE)
 	  vec <- prcomp(xtrain, center = TRUE, scale = TRUE)$rotation
      z <- xtrain %*% vec  ## PCA scores
      xnew <- t( ( t(xtest) - mx ) / s ) ## standardize the xnew values

      for ( j in 1:maxk ) {
        zzk <- crossprod(z[, 1:j])
        be <- vec[, 1:j] %*% solve( zzk, crossprod( z[, 1:j], ytrain ) )
        ## b is the PCA based coefficients
        est <- as.vector( m + xnew %*% be )  ## predicted values for PCA model
        er[j] <- sum( (ytest - est)^2 ) / rmat
      }
      return(er)
    }
    stopCluster(cl)

    runtime <- proc.time() - runtime
  }

  mspe <- Rfast::colmeans(msp)
  bias <- msp[, which.min(mspe)] - Rfast::rowMins(msp, value = TRUE)   ## apply(msp, 1, min)  ## TT estimate of bias
  estb <- mean( bias )  ## TT estimate of bias

  if ( graph )  plot(1:maxk, mspe, xlab = "Number of principal components", ylab = "MSPE", type = "b")
  names(mspe) <- paste("PC", 1:maxk, sep = " ")
  performance <- c( min(mspe) + estb, estb)
  names(performance) <- c("MSPE", "Estimated bias")
  list(msp = msp, mspe = mspe, k = which.min(mspe), performance = performance, runtime = runtime)
}
