pprcomp.tune <- function(y, x, nfolds = 10, folds = NULL, seed = FALSE, nterms = 1:10, type = "alr") {

  runtime <- proc.time()
  if ( type == "alr" ) {
    x <- alr(x)
  } else  x <- Rfast::Log(x)
  x <- as.data.frame(x)

  n <- dim(x)[1]
  lnt <- length(nterms)
  ina <- 1:n
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                                           stratified = FALSE, seed = seed)
  nfolds <- length(folds)
  names <- paste("fold", 1:nfolds)
  mse <- matrix(nrow = nfolds, ncol = lnt)
  rownames(mse) <- paste("Fold", 1:nfolds, sep = " ")
  colnames(mse) <- paste("nterms=", nterms, sep = " ")

  for (i in 1:nfolds) {
    nu <- folds[[ i ]]
    ytrain <- y[-nu]
    xtrain <- x[-nu, ]
    ytest <- y[nu]
    xtest <- x[nu, ]
    for (j in 1:lnt ) {
      mod <- ppr( ytrain ~., data = xtrain, nterms = nterms[j] )
      est <- predict(mod, newdata = xtest)
      mse[i, j] <- mean( (ytest - est)^2 )
    }
  }

  runtime <- proc.time() - runtime
  list(runtime = runtime, mse = mse )
}
