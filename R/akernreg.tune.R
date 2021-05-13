## Response is compositional data
akernreg.tune <- function(y, x, a = seq(0.1, 1, by = 0.1), h = seq(0.1, 1, length = 10),
                          nfolds = 10, folds = NULL, seed = FALSE) {
  runtime <- proc.time()
  D <- dim(y)[2]
  n <- dim(y)[1]
  if ( min(y) == 0 )  a <- a[a > 0]
  la <- length(a)
  ina <- 1:n
  if ( !is.matrix(x) )  x <- as.matrix(x)
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds, stratified = FALSE, seed = seed)
  nfolds <- length(folds)
  names <- paste("fold", 1:nfolds)
  est <- sapply(names, function(x) NULL)
  kl <- list()
  nh <- length(h)
  kl[[ 1 ]] <- matrix(NA, nrow = la, ncol = nh)
  rownames( kl[[ 1 ]] ) <- paste("alpha", a, sep = " ")
  colnames( kl[[ 1 ]] ) <- paste("h=", h, sep = "")
  js <- kl

  for (i in 1:nfolds) {
    nu <- folds[[ i ]]
    ytrain <- y[-nu, , drop = FALSE]
    xtrain <- x[-nu, , drop = FALSE]
    ytest <- y[nu, , drop = FALSE]
    xtest <- x[nu, , drop = FALSE]
    est[[ i ]] <- akern.reg(xnew = xtest, y = ytrain, x = xtrain, a = a, h = h)
    kl[[ i ]] <- kl[[ 1 ]]
    js[[ i ]] <- js[[ 1 ]]
    for (j in 1:la)  {
      for (vim in 1:nh) {
        zest <- est[[ i ]][[ j ]][[ vim ]]
        ela <- abs( ytest * log( ytest / zest ) )
        ela[ is.infinite(ela) ] <- NA
        kl[[ i ]][j, vim] <- 2 * mean(ela , na.rm = TRUE)
        ela2 <- ytest * log( 2 * ytest / (ytest + zest) ) + zest * log( 2 * zest / (ytest + zest) )
        ela2[ is.infinite(ela2) ] <- NA
        js[[ i ]][j, vim] <- mean(ela2, na.rm = TRUE)
      }
    }
  }

  names(kl) <- names(js) <- paste("Folds", 1:nfolds, sep = " ")
  jsula <- kula <- 0
  for (i in 1:nfolds)  {
    kula <- kula + kl[[ i ]]
    jsula <- jsula + js[[ i ]]
  }
  kula <- kula / nfolds
  jsula <- jsula / nfolds
  bkl <- which( kula == min(kula, na.rm = TRUE), arr.ind = TRUE )
  bjs <- which( jsula == min(jsula, na.rm = TRUE), arr.ind = TRUE )

  res <- list( kl = kula, js = jsula, klmin = min(kula, na.rm = TRUE), jsmin = min(jsula, na.rm = TRUE),
               kl.alpha = a[ bkl[1] ], kl.h = h[ bkl[2] ], kl.alpha = a[ bjs[1] ], kl.h = h[ bjs[2] ] )

  runtime <- proc.time() - runtime
  res$runtime <- runtime
  res
}



