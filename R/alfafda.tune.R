################################
#### Classification for compositional data using the alpha-transformation
#### Tuning the hyper-parameters via K-fold cross-validation
#### Tsagris Michail 7/2015
#### References: Tsagris, M., Preston S. and Wood A.T.A. (2016).
#### Improved classication for compositional data using the alpha-transformation
#### Journal of Classification (To appear)
#### http://arxiv.org/pdf/1506.04976v2.pdf
#### mtsagris@yahoo.gr
################################
alfafda.tune <- function(x, ina, a = seq(-1, 1, by = 0.1), nfolds = 10, folds = NULL,
                         stratified = TRUE, seed = FALSE, graph = FALSE) {
  ## x contains the compositonal data
  ## ina is the grouping variable
  ## a is the grid of values of a
  ## M is th number of folds
  ## ncores is the number of cores to be used
  ## if mat is NULL the folds happen internally
  ## if you already have folds, provide the indices of the data
  ## in a matrix form, each column corresponds to a fold
  if ( min(x) == 0 )  a = a[ a > 0 ]   ## if you have zero values, only positive alphas are allowed
  info <- list()
  ina <- as.numeric(ina)
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                                           stratified = stratified, seed = seed)
  per <- matrix(nrow = length(a), ncol = nfolds)
  runtime <- proc.time()
  for ( i in 1:length(a) ) {
    z <- Compositional::alfa(x, a[i])$aff  ## apply the alpha-transformation
    for (k in 1:nfolds) {
      test <- z[ folds[[ k ]], , drop = FALSE ]   ## test sample
      id <- ina[ folds[[ k ]] ] ## groups of test sample
      train <- z[ -folds[[ k ]], , drop = FALSE]  ## training sample
      ida <- ina[ -folds[[ k ]] ]   ## groups of training sample
      mod <- mda::fda(ida ~ train)
      g <- predict(mod, test)
      per[i, k] <- mean( g == id )
    }
  }
  runtime <- proc.time() - runtime
  performance <- Rfast::rowmeans(per)
  names(performance) <- paste("alfa=", a, sep = "")

  if ( graph ) plot(a, performance, type = "l", ylim = c( min(performance), max(performance) ),
                    ylab = "Estimated performance", xlab = expression(paste(alpha, " values")), cex.lab = 1.3 )

  list( per = per, performance = performance, best_a = a[ which.max(per) ], runtime = runtime )
}

