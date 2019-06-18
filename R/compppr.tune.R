compppr.tune <- function(y, x, nfolds = 10, folds = NULL, seed = FALSE, nterms = 1:10, type = "alr", 
  yb = NULL, B = 1000 ) {
  
  runtime <- proc.time() 
  if ( is.null(yb) )  {
    if ( type == "alr" ) {
      yb <- alr(y)
    } else   yb <- alfa(y, 0, h = TRUE)$aff
  }
  
  D <- dim(y)[2]
  n <- dim(y)[1]
  nt <- nterms
  if (seed)  set.seed(12345)
  ind <- sample(n, n)
  if ( is.null(folds) ) {
    nf <- round( n / nfolds ) 
    folds <- list()
    for ( i in 1:c(nfolds - 1) )  {
      folds[[ i ]] <- ind[1:nf]  
      ind <- ind[ - c(1:nf) ]
    }
    folds[[ nfolds ]] <- ind      
  } 

  names <- paste("fold", 1:nfolds)
  est <- sapply(names, function(x) NULL)
  kl <- matrix( NA, nrow = nfolds, ncol = length(nt) )

  for (i in 1:nfolds) {
    nu <- folds[[ i ]]
    ytrain <- y[-nu, ]
    ybtrain <- yb[-nu, ]
    xtrain <- x[-nu, ]
    ytest <- y[nu, ]
    ybtest <- yb[nu, ]
    xtest <- x[nu, ]
    for (j in 1:length(nt) ) {
      mod <- Compositional::comp.ppr(ytrain, xtrain, nterms = nt[j], type = "alr", xnew = xtest, yb = ybtrain)
      est[[ i ]][[ j ]] <- mod$est
      kl[i, j] <- mean( abs( ytest * log(ytest/mod$est) ), na.rm = TRUE)
    }
  }
  rownames(kl) <- paste("Fold", 1:nfolds, sep = " ")
  colnames(kl) <- paste("nterms=", nt, sep = " ")

  for ( j in nt ) {
    for ( i in 2:nfolds ) {
      est[[ i ]][[ j ]] <- rbind( est[[ i - 1 ]][[ j ]], est[[ i ]][[ j ]] )
    }
  }

  for (j in nt)   est[[ j ]] <- est[[ nfolds ]][[ j ]] ## compositional predictions
  if ( max(nt) < nfolds )  for (j in nfolds:c(max(nt) + 1) )  est[[ j ]] <- NULL
  ina <- unlist(folds)
  yy <- y[ina, ]
  
  if ( B > 1 ) {
    perf <- numeric(B)
    for (i in 1:B) {
      ind <- sample(n, n, replace = TRUE)
      a <- numeric( length(nt) )
      for ( j in 1:length(nt) )  a[j] <- mean( abs( yy[ind, ] * log(yy[ind, ]/est[[ j ]][ind, ]) ), na.rm = TRUE)
      m <- which.min(a)
      perf[i] <- mean( abs( yy[-ind, ] * log(yy[-ind, ]/est[[ m ]][-ind, ]) ), na.rm = TRUE)
    }
	res <- list( kl = 2 * kl, bc.perf = 2 * mean(perf) )
  }  else 	res <- list( kl = 2 * kl, bc.perf = 2 * kl )
  
  runtime <- proc.time() - runtime
  res$runtime <- runtime
  res  
}
