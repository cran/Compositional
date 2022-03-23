compppr.tune <- function(y, x, nfolds = 10, folds = NULL, seed = NULL, nterms = 1:10, type = "alr", yb = NULL) {

  runtime <- proc.time()
  if ( is.null(yb) )  {
    if ( type == "alr" ) {
      yb <- alr(y)
    } else   yb <- alfa(y, 0, h = TRUE)$aff
  }

  n <- dim(y)[1]
  D <- dim(y)[2]
  nt <- nterms
  ina <- 1:n
  if ( !is.matrix(x) )  x <- as.matrix(x)
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                  stratified = FALSE, seed = seed)
  nfolds <- length(folds)

  names <- paste("fold", 1:nfolds)
  est <- sapply(names, function(x) NULL)
  kl <- matrix( NA, nrow = nfolds, ncol = length(nt) )
  rownames(kl) <- paste("Fold", 1:nfolds, sep = " ")
  colnames(kl) <- paste("nterms=", nt, sep = " ")

  for (i in 1:nfolds) {
    nu <- folds[[ i ]]
    ytrain <- y[-nu, ]
    ybtrain <- yb[-nu, ]
    xtrain <- x[-nu, , drop = FALSE]
    ytest <- y[nu, ]
    ybtest <- yb[nu, ]
    xtest <- x[nu, , drop = FALSE]
    for (j in 1:length(nt) ) {
      mod <- Compositional::comp.ppr(ytrain, xtrain, nterms = nt[j], type = "alr", xnew = xtest, yb = ybtrain)
      est[[ i ]][[ j ]] <- mod$est
      kl[i, j] <- mean( abs( ytest * log(ytest/mod$est) ), na.rm = TRUE)
    }
  }

 	res <- list( kl = 2 * kl )

  runtime <- proc.time() - runtime
  res$runtime <- runtime
  res
}
