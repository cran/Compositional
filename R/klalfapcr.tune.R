klalfapcr.tune <- function(y, x, covar = NULL, M = 10, maxk = 50,
                           a = seq(-1, 1, by = 0.1), mat = NULL, graph = FALSE,
                           tol = 1e-07, maxiters = 50) {
  n <- dim(x)[1]
  p <- dim(x)[2] - 1
  if ( min(x) == 0 )  a <- a[ a > 0 ]
  if ( maxk > p )   maxk <- p
  if ( !is.null(covar) )  covar <- as.matrix(covar)
  if ( is.null(mat) ) {
    nu <- sample(1:n, min( n, round(n / M) * M ) )
    options(warn = -1)
    mat <- matrix( nu, ncol = M )
  } else  mat <- mat
  M <- ncol(mat)
  mspe <- list()
  msp <- matrix( nrow = M, ncol = maxk )
  colnames(msp) <- paste("PC", 1:maxk, sep = " ")
  for ( i in 1:length(a) ) {
    xa <- Compositional::alfa(x, a[i])$aff
    for (vim in 1:M) {
      ytest <- y[ mat[, vim], , drop = FALSE ]
      ytrain <- y[ -mat[, vim], , drop = FALSE ]
      xtrain <- xa[ -mat[, vim],,  drop = FALSE ]
      xtest <- xa[ mat[, vim], , drop = FALSE ]
      com <- sum(ytest * log(ytest), na.rm = TRUE)
      mod <- prcomp(xtrain, center = FALSE)
      vec <- mod$rotation
      za <- mod$x
      zanew <- xtest %*% vec
      for ( j in 1:maxk ) {
        if ( !is.null(covar) ) {
          z <- cbind(za[, 1:j, drop = FALSE],
                     covar[ -mat[, vim], drop = FALSE ] )
          znew <- cbind(zanew[, 1:j, drop = FALSE],
                        covar[ mat[, vim], drop = FALSE ] )
        } else {
          z <- za[, 1:j, drop = FALSE ]
          znew <- zanew[, 1:j, drop = FALSE]
        }
        est <- Compositional::kl.compreg(y = ytrain, x = z, xnew = znew,
                                         tol = 1e-07, maxiters = maxiters)$est
        res <- sum(ytest * log(est), na.rm = TRUE)
        msp[vim, j] <- com - res * is.finite(res)
      }
    }
    mspe[[ i ]] <- msp
  }
  names(mspe) <- paste("alpha=", a, sep = "")
  performance <- lapply(mspe, colMeans)
  performance <- matrix( unlist(performance), ncol = maxk, byrow = TRUE )
  colnames(performance) <- paste("PC", 1:maxk, sep = " ")
  rownames(performance) <- paste("alpha", a, sep = " ")
  poia <- which(performance == min(performance, na.rm = TRUE), arr.ind = TRUE)
  params <- c( a[ poia[, 1] ], poia[, 2] )
  names(params) <- c("best alpha", "best k")
  if ( graph ) {
    filled.contour(a, 1:maxk, performance, xlab = expression(paste(alpha, " values")),
                   ylab = "Number of PCs", cex.lab = 1.3 )
  }
  list(mspe = mspe, performance = performance, best.perf = min(performance), params = params)
}
