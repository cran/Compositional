cv.dda <- function(x, ina, nfolds = 10, folds = NULL, stratified = TRUE, seed= FALSE) {
  runtime <- proc.time()
  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds, stratified = stratified, seed = seed)
  nfolds <- length(folds)
  per <- numeric(nfolds)

  for ( k in 1:nfolds ) {
    xtrain <- x[ -folds[[ k ]], , drop = FALSE]  ## training sample
    idtrain <- ina[ -folds[[ k ]] ]   ## groups of training sample
    xtest <- x[ folds[[ k ]], , drop = FALSE ]   ## test sample
    idtest <- ina[ folds[[ k ]] ] ## groups of test sample
    mat <- matrix(nrow = length(idtest), ncol = g)
    for (j in 1:g) {
      a <- Compositional::diri.nr(xtrain[idtrain == j, ])$param
      mat[, j] <- Compositional::ddiri(xtest, a, logged = TRUE)
    }
    est <- Rfast::rowMaxs(mat)
    per[k] <- mean(idtest == est)
  }  ##  end  for (k in 1:nfolds)

  percent <- mean(per)
  runtime <- proc.time() - runtime
  list(percent = percent, runtime = runtime)
}
