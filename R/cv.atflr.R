cv.atflr <- function(y, x, a = seq(0.1, 1, by = 0.1), nfolds = 10, folds = NULL, seed = NULL) {

  if ( min(y) == 0 )  a <- abs(a)
  n <- dim(y)[1]
  ina <- 1:n
  if (is.null(folds))  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                                         stratified = FALSE, seed = seed)
  nfolds <- length(folds)
  la <- length(a)
  js <- kl <- matrix(nrow = nfolds, ncol = la)

  runtime <- proc.time()
  for ( j in 1:la ) {
    ya <- y^a[j]
    ya <- ya / Rfast::rowsums(ya)
    for ( i in 1:nfolds ) {
      ytest <- y[ folds[[ i ]], ]  ## test set dependent vars
      ytrain <- ya[ -folds[[ i ]], ]  ## train set dependent vars
      xtest <- x[ folds[[ i ]], ]  ## test set independent vars
      xtrain <- x[ -folds[[ i ]], ]  ## train set independent vars
      est <- Compositional::tflr(ytrain, xtrain, xnew = xtest)$est^( 1 / a[j] )
      est <- est / Rfast::rowsums(est)
      ela <- abs( ytest * log( ytest / est ) )
      ela[ is.infinite(ela) ] <- NA
      kl[i, j] <-  2 * mean(ela, na.rm = TRUE)
      ela2 <- ytest * log( 2 * ytest / (ytest + est) ) + est * log( 2 * est / (ytest + est) )
      ela2[ is.infinite(ela2) ] <- NA
      js[i, j] <- mean(ela2, na.rm = TRUE)
    }
  }
  runtime <- proc.time() - runtime

  kl <- Rfast::colmeans(kl)
  js <- Rfast::colmeans(js)
  names(kl) <- names(js) <- paste("alpha=", a, sep = "")
  list(runtime = runtime, kl = kl, js = js )
}
