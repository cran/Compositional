cv.olscompcomp <- function(y, x, tol = 1e-4, nfolds = 10, folds = NULL, seed = NULL) {

  n <- dim(y)[1]
  ina <- 1:n
  if (is.null(folds))  folds <- Compositional::makefolds(ina, nfolds = nfolds, stratified = FALSE, seed = seed)
  nfolds <- length(folds)
  js <- kl <- numeric(nfolds)

  runtime <- proc.time()
  for ( i in 1:nfolds) {
    ytest <- y[ folds[[ i ]], ]  ## test set dependent vars
    ytrain <- y[ -folds[[ i ]], ]  ## train set dependent vars
    xtest <- x[ folds[[ i ]], -1, drop = FALSE]  ## test set independent vars
    xtrain <- x[ -folds[[ i ]], -1, drop = FALSE]  ## train set independent vars
    est <- Compositional::ols.compcomp(ytrain, xtrain, xnew = xtest)$est
    ela <- abs( ytest * log( ytest / est ) )
    ela[ is.infinite(ela) ] <- NA
    kl[i] <-  2 * mean(ela, na.rm = TRUE)
    ela2 <- ytest * log( 2 * ytest / (ytest + est) ) + est * log( 2 * est / (ytest + est) )
    ela2[ is.infinite(ela2) ] <- NA
    js[i] <- mean(ela2, na.rm = TRUE)
  }
  runtime <- proc.time() - runtime

  perf <- c( mean(kl), mean(js) )
  names(perf) <- c( "KL", "JS")
  list(runtime = runtime, kl = kl, js = js, perf = perf )
}
