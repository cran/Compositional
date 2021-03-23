## Predictors are compositional data
alfaknnreg.tune <- function(y, x, a = seq(-1, 1, by = 0.1), k = 2:10, nfolds = 10,
apostasi = "euclidean", method = "average", folds = NULL, seed = FALSE, graph = FALSE) {
  if ( min(x) == 0 )  a <- a[a>0]  ## checks for any zeros in the data
  n <- dim(x)[1]  ## sample size
  ina <- 1:n
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                                           stratified = FALSE, seed = seed)
  nfolds <- length(folds)
  runtime <- proc.time()
  mspe <- matrix(nrow = length(a), ncol = length(k) )
  colnames(mspe) <- paste("k=", k, sep = "")
  rownames(mspe) <- paste("alpha=", a, sep = "")
  for (i in 1:length(a) ) {
    z <- alfa(x, a[i], h = FALSE)$aff
    mspe[i, ] <- Rfast::knn.cv(folds = folds, nfolds = nfolds, y = y, x = z, k = k,
                 dist.type = apostasi, type = "R", method = method)$crit
  }
  runtime <- proc.time() - runtime
  if ( graph )  filled.contour(a, k, mspe, col = grey(1:11/11), ylab = "k nearest-neighbours",
                       xlab = expression(paste(alpha, " values")), cex.lab = 1.2, cex.axis = 1.2)
  opt <- min(mspe)
  confa <- as.vector( which(mspe == opt, arr.ind = TRUE)[1, ] )
  list( mspe = mspe, performance = min(mspe), opt_a = a[ confa[1] ],
        opt_k = confa[2] + 1, runtime = runtime )
  }

