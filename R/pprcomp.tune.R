pprcomp.tune <- function(y, x, nfolds = 10, folds = NULL, seed = NULL,
                         nterms = 1:10, type = "log", graph = FALSE) {

  runtime <- proc.time()
  if ( type == "alr" ) {
    x <- Compositional::alr(x)
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
  mse <- Rfast::colmeans(mse)

  if ( graph ) {
    plot(nterms, mse, type = "b", xlab = "Number of terms", pch = 19, cex.lab = 1.3, cex.axis = 1.3,
      ylab = "Mean squared error of prediction", lwd = 2, col = "green")
    abline(v = nterms, lty = 2, col = "lightgrey")
    abline(h = seq(min(mse), max(mse), length = 10), lty = 2, col = "lightgrey" )
  }

  list(runtime = runtime, mse = mse, opt.nterms = nterms[ which.min(mse) ], performance = min(mse) )
}


