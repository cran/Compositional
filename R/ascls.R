ascls <- function(y, x, a = seq(0.1, 1, by = 0.1), xnew ) {

  if ( min(y) == 0 )  a <- abs(a)
  py <- dim(y)[2]   ;    px <- dim(x)[2]
  pyx <- py * px    ;    n <- dim(y)[1]

  xx <- crossprod(x)
  XX <- matrix(0, pyx, pyx)
  ind <- matrix( 1:pyx, ncol = px, byrow = TRUE)
  for ( i in 1:py )  XX[ ind[i, ], ind[i, ] ] <- xx
  A <- matrix(0, pyx, pyx)
  for ( i in 1:px )  A[i, ind[, i]] <- 1
  A <- t( rbind( A, diag(pyx), -diag(pyx) ) )
  A <- A[, -c( (px + 1): pyx) ]
  bvec <- c( rep(1, px), rep(0, pyx), rep(-1, pyx) )

  res <- list()
  for ( j in 1:length(a) ) {
    ya <- y^a[j]
    ya <- ya / Rfast::rowsums(ya)
    dvec <- as.vector( crossprod(x, ya) )
    f <- try( quadprog::solve.QP( Dmat = XX, dvec = dvec, Amat = A, bvec = bvec,
                                  meq = px ), silent = TRUE )
    if ( identical(class(f), "try-error") ) {
      f <- quadprog::solve.QP( Dmat = Matrix::nearPD(XX)$mat, dvec = dvec, Amat = A, bvec = bvec, meq = px )
    }
    be <- matrix( abs(f$solution), ncol = py)
    est <- ( xnew %*% be )^( 1/a[j] )
    est <- est / Rfast::rowsums(est)
    res[[ j ]] <- est
  }

  res
}

