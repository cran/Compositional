################################
#### Principal components regression for binary and poisson regression
#### Selection of the number of principal components
#### via K-fold cross validation
#### Tsagris Michail 1/2016
#### mtsagris@yahoo.gr
################################

glmpcr.tune <- function(y, x, M = 10, maxk = 10, oiko = "binomial", seed = FALSE,
                        ncores = 2, graph = TRUE) {
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
  x <- as.matrix(x)
  y <- as.vector(y)
  n <- nrow(x)
  p <- ncol(x)
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
        est <- glm.pcr(ytrain, xtrain, j, oiko = oiko, xtest)$est
        if (oiko == "binomial") {
          ri <-  -2 *( ytest * log(est) + (1 - ytest) * log(1 - est) )
        } else {
          ri <- 2 * ( ytest * log(ytest / est) )
        }
        msp[vim, j] <- sum( ri )
      }
    }
  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    er <- numeric(maxk)
    msp <- foreach(vim = 1:M, .combine = rbind, .export = "glm.pcr") %dopar% {
      ## will always be the same
      ytest <- as.vector( y[mat[, vim] ] )  ## test set dependent vars
      ytrain <- as.vector( y[-mat[, vim] ] )  ## train set dependent vars
      xtrain <- as.matrix( x[-mat[, vim], ] )  ## train set independent vars
      xtest <- as.matrix( x[mat[, vim], ] )  ## test set independent vars
      for ( j in 1:maxk ) {
        est <- glm.pcr(ytrain, xtrain, j, oiko = oiko, xtest)$est
        if (oiko == "binomial") {
          ri <-  -2 * ( ytest * log(est) + (1 - ytest) * log(1 - est) )
        } else {
          ri <- 2 * ( ytest * log( ytest / est ) )
        }
        er[j] <- sum( ri )
      }
      return(er)
    }
    stopCluster(cl)
  }
  mpd <- colMeans(msp)
  bias <- msp[ ,which.min(mpd)] - apply(msp, 1, min)  ## TT estimate of bias
  estb <- mean( bias )  ## TT estimate of bias
  if (graph == TRUE) {
    plot(1:maxk, mpd, xlab = "Number of principal components",
    ylab = "Mean predicted deviance", type = "b")
  }
  names(mpd) <- paste("PC", 1:maxk, sep = " ")
  performance <- c( min(mpd) + estb, estb)
  names(performance) <- c("MPD", "Estimated bias")
  list(msp = msp, mpd = mpd, k = which.min(mpd), performance = performance)
}
