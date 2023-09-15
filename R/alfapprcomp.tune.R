alfapprcomp.tune <- function(y, x, nfolds = 10, folds = NULL, seed = NULL,
                         nterms = 1:10, a = seq(-1, 1, by = 0.1), graph = FALSE) {

  runtime <- proc.time()

  if ( min(x) == 0 )  a <- a[ a > 0 ]
  n <- dim(x)[1]
  lnt <- length(nterms)
  ina <- 1:n
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                                           stratified = FALSE, seed = seed)
  nfolds <- length(folds)
  names <- paste("fold", 1:nfolds)
  intmse <- matrix(nrow = nfolds, ncol = lnt)
  mse <- matrix(nrow = length(a), ncol = lnt)
  rownames(mse) <- paste("alpha=", a, sep = " ")
  colnames(mse) <- paste("nterms=", nterms, sep = " ")

  for ( k in 1:length(a) ) {
    z <- Compositional::alfa(x, a[k])$aff  ## apply the alpha-transformation
    z <- as.data.frame(z)
    for (i in 1:nfolds) {
      nu <- folds[[ i ]]
      ytrain <- y[-nu]
      ztrain <- z[-nu, ]
      ytest <- y[nu]
      ztest <- z[nu, ]
      for (j in 1:lnt ) {
        mod <- ppr( ytrain ~., data = ztrain, nterms = nterms[j] )
        est <- predict(mod, newdata = ztest)
        intmse[i, j] <- mean( (ytest - est)^2 )
      }
    }
    mse[k, ] <- Rfast::colmeans(intmse)
  }

  runtime <- proc.time() - runtime

  if ( graph ) {
    if ( graph )  filled.contour( a, nterms, ylab = "Number of terms", cex.lab = 1.2,
                                  cex.axis = 1.2, xlab = expression(paste(alpha, " values") ) )
  }

  pou <- which( mse == min(mse), arr.ind = TRUE )

  list(runtime = runtime, mse = mse, opt.nterms = nterms[ pou[2] ], opt.alpha = a[ pou[1] ], performance = min(mse) )
}


