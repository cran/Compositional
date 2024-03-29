###############################
#### Classification for compositional data using the alpha-transformation
#### Tuning the k-NN algorithm
#### Tsagris Michail 8/2015
#### References: Tsagris, M., Preston S. and Wood A.T.A. (2016).
#### Improved classication for compositional data using the alpha-transformation
#### Journal of Classification (To appear)
#### http://arxiv.org/pdf/1506.04976v2.pdf
#### mtsagris@yahoo.gr
################################
alfaknn.tune <- function(x, ina, nfolds = 10, k = 2:5, mesos = TRUE, a = seq(-1, 1, by = 0.1),
                         apostasi = "euclidean", rann = FALSE, folds = NULL, stratified = TRUE, seed = NULL, graph = FALSE) {
  if ( min(x) == 0 )  a <- a[a > 0]  ## checks for any zeros in the data
  n <- dim(x)[1]  ## sample size
  ina <- as.numeric(ina) ## makes sure ina is numeric
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                                           stratified = stratified, seed = seed)
  nfolds <- length(folds)
  runtime <- proc.time()
  if ( !rann ) {
    ela <- matrix( nrow = length(a), ncol = length(k) )
    for (i in 1:length(a) ) {
      z <- Compositional::alfa(x, a[i], h = FALSE)$aff
      ela[i, ] <- Rfast::knn.cv(folds = folds, nfolds = nfolds, y = ina, x = z, k = k, dist.type = apostasi, type = "C", freq.option = 1)$crit
    }

  } else {
    per <- array( dim = c( nfolds, length(k), length(a) ) )  ## The estimated percentages
    for ( i in 1:length(a) ) {
      z <- Compositional::alfa(x, a[i], h = FALSE)$aff
      for (vim in 1:nfolds) {
        id <- ina[ folds[[ vim ]] ]   ## groups of test sample
        ina2 <- ina[ -folds[[ vim ]] ]   ## groups of training sample
        aba <- as.vector( folds[[ vim ]] )
        aba <- aba[aba > 0]
        g <- Compositional::alfa.knn(z[aba, ], z[-aba, ], ina = ina2, a = NULL, k = k, rann = TRUE)
        per[vim, , i] <- Rfast::colmeans( g == id )
      }
    }
    ela <- t( colMeans(per) )
  }
  colnames(ela) <- paste("k=", k, sep = "")
  rownames(ela) <- paste("alpha=", a, sep = "")
  runtime <- proc.time() - runtime
  if ( graph )  filled.contour( a, k, ela, ylab = "k nearest-neighbours", cex.lab = 1.2, cex.axis = 1.2,
                                xlab = expression(paste(alpha, " values") ) )
  opt <- max(ela)
  confa <- as.vector( which(ela == opt, arr.ind = TRUE)[1, ] )
  list( ela = ela, performance = max(ela), best_a = a[ confa[1] ], best_k = confa[2] + 1, runtime = runtime )

}
