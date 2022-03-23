cv.compnb <- function(x, ina, type = "beta", folds = NULL, nfolds = 10,
             stratified = TRUE, seed = NULL, pred.ret = FALSE) {

  ina <- as.numeric(ina)
  if ( is.null(folds) ) {
    folds <- Compositional::makefolds(ina, nfolds = nfolds, stratified = stratified, seed = seed)
  }
  nfolds <- length(folds)
  crit <- numeric(nfolds)
  preds <- NULL
  if (pred.ret) {
    names <- paste("Fold", 1:nfolds)
    preds <- sapply(names, function(x) NULL)
  }

  for (i in 1:nfolds) {
    inatrain <- ina[ -folds[[ i ]] ]
    xtrain <- x[-folds[[ i ]], ]
    inatest <- ina[ folds[[ i ]] ]
    xtest <- x[folds[[ i ]], ]
    est <- Compositional::comp.nb(xnew = xtest, x = xtrain, ina = inatrain, type = type)$est
    if ( pred.ret )   preds[[i]] <- est
    crit[i] <- mean(est == inatest)
  }

  list(preds = preds, crit = mean(crit))

}
