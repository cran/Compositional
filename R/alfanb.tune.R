alfanb.tune <- function(x, ina, a = seq(-1, 1, by = 0.1), type = "gaussian",
               folds = NULL, nfolds = 10, stratified = TRUE, seed = FALSE) {

  ina <- as.numeric(ina)
  if ( is.null(folds) ) {
    folds <- Compositional::makefolds(ina, nfolds = nfolds, stratified = stratified, seed = seed)
  }
  nfolds <- length(folds)
  la <- length(a)
  mat <- numeric(la)
  names(mat) <- paste("a=", a, sep = "")
  if (la > 1) {
    if ( min(x) == 0 )  a <- a[a > 0]
  }

  if ( type == "gaussian" ) {
    nb <- Rfast::gaussian.nb
  } else if ( type == "cauchy" ) {
    nb <- Rfast2::cauchy.nb
  } else if ( type == "laplace" ) {
    nb <- Rfast2::laplace.nb
  }

  for (j in 1:la) {
    y <- Compositional::alfa(x, a[j])$aff
    crit <- numeric(nfolds)

    for (i in 1:nfolds) {
      inatrain <- ina[ -folds[[ i ]] ]
      ytrain <- y[-folds[[ i ]], ]
      inatest <- ina[ folds[[ i ]] ]
      ytest <- y[folds[[ i ]], ]
      est <- nb(xnew = ytest, x = ytrain, ina = inatrain)$est
      crit[i] <- mean(est == inatest)
    }
    mat[j] <- mean(crit)
  }

  mat
}
