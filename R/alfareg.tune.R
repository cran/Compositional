################################
#### Tuning the alfa in alfa-regression via K-fold cross validation
#### Tibshirani and Tibshirani method
#### Tsagris Michail 11/2015
#### mtsagris@yahoo.gr
#### References: Tsagris Michail (2015)
#### Regression analysis with compositional data containing zero values
#### Chilean Journal of Statistics, 6(2): 47-57
################################
alfareg.tune <- function(y, x, a = seq(0.1, 1, by = 0.1), nfolds = 10, folds = NULL, nc = 1,
                         seed = NULL, graph = FALSE) {
  ## y is the compositional data (dependent variable)
  ## x is the independent variables
  ## a is a range of values of alpha
  ## nfolds is the number of folds for the K-fold cross validation
  ## nc is how many cores you want to use, default value is 2
  if ( min(y) == 0 )  a <- a[a>0]
  la <- length(a)
  n <- dim(y)[1]
  ina <- 1:n
  x <- model.matrix(y ~., data.frame(x) )
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                                           stratified = FALSE, seed = seed)
  nfolds <- length(folds)

  if (nc <= 1) {
    apa <- proc.time()
    kula <- matrix(nrow = nfolds, ncol = la)
    for (j in 1:la) {
      ytr <- Compositional::alfa(y, a[j])$aff
      for (i in 1:nfolds) {
        xu <- x[ folds[[ i ]], -1 , drop = FALSE]
        yu <- y[ folds[[ i ]], ]
        xa <- x[ -folds[[ i ]], -1]
        yb <- ytr[ -folds[[ i ]], ]
        mod <- Compositional::alfa.reg(yu, xa, a[j], xnew = xu, yb = yb)
        yest <- mod$est
        kula[i, j] <- 2 * mean(yu * log(yu / yest), na.rm = TRUE)
      }
    }

    kl <- Rfast::colmeans(kula)
    opt <- a[ which.min(kl) ]
    val <- which.min(kl)
    per <- min(kl, na.rm = TRUE)
    pera <- Rfast::rowMins(kula, value = TRUE)  ## apply(kula, 1, min)
    apa <- proc.time() - apa

  } else {
    apa <- proc.time()
    val <- matrix(a, ncol = nc) ## if the length of a is not equal to the
    ## dimensions of the matrix val a warning message should appear
    cl <- parallel::makePSOCKcluster(nc)
    doParallel::registerDoParallel(cl)
    if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                                             stratified = FALSE, seed = seed)
    kula <- foreach::foreach(j = 1:nc, .combine = cbind, .packages = "Rfast", .export = c("alfa.reg",
	        "alfa", "helm", "comp.reg", "multivreg", "rowsums", "colmeans", "colVars") ) %dopar% {
       ba <- val[, j]
       ww <- matrix(nrow = nfolds, ncol = length(ba) )
       for ( l in 1:length(ba) ) {
          ytr <- Compositional::alfa(y, ba[l])$aff
          for (i in 1:nfolds) {
            xu <- x[ folds[[ i ]], -1 , drop = FALSE]
            yu <- y[ folds[[ i ]], ]
            xa <- x[ -folds[[ i ]], -1]
            yb <- ytr[ -folds[[ i ]], ]
            mod <- alfa.reg(yu, xa, ba[l], xnew = xu, yb = yb)
            yest <- mod$est
            ww[i, l] <- 2 * mean(yu * log(yu / yest), na.rm = TRUE)
          }
       }
       return(ww)
    }
    parallel::stopCluster(cl)

    kula <- kula[, 1:la]
    kl <- Rfast::colmeans(kula)
    opt <- a[ which.min(kl) ]
    val <- which.min(kl)
    per <- min(kl, na.rm = TRUE)
    pera <- Rfast::rowMins(kula, value = TRUE)   ## apply(kula, 1, min)
    apa <- proc.time() - apa
  }

  if ( graph ) {
    plot( a, kl, type = 'b', ylim = c( min(kl), max(kl) ), xlab = expression(alpha),
	      ylab = '2 * Kullback Leibler divergence', cex.lab = 1.2, cex.axis = 1.2, pch = 16,
	      col = "green", lwd = 2 )
    abline(v = a, col = "lightgrey", lty = 2)
    abline(h = seq(min(kl), max(kl), length = 10), col = "lightgrey", lty = 2)
  }

  list(runtime = apa, kula = kula, kl = kl, opt = opt, value = per)
}
