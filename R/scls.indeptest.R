scls.indeptest <- function(y, x, R = 999) {

  py <- dim(y)[2]   ;    px <- dim(x)[2]
  pyx <- py * px    ;    n <- dim(y)[1]

  dvec <- as.vector( crossprod(x, y) )
  xx <- crossprod(x)
  Dmat <- matrix(0, pyx, pyx)
  ind <- matrix( 1:pyx, ncol = px, byrow = TRUE )
  for ( i in 1:py )  Dmat[ ind[i, ], ind[i, ] ] <- xx
  A <- matrix(0, pyx, pyx)
  for ( i in 1:px )  A[i, ind[, i]] <- 1
  A <- t( rbind( A, diag(pyx), -diag(pyx) ) )
  A <- A[, -c( (px + 1): pyx) ]
  bvec <- c( rep(1, px), rep(0, pyx), rep(-1, pyx) )

  if ( det(Dmat) < 1e-20 )  Dmat <- Matrix::nearPD(Dmat)$mat

  f <- quadprog::solve.QP( Dmat = Dmat, dvec = dvec, Amat = A, bvec = bvec, meq = px )
  mse <- f$value
  pmse <- numeric(R)
  for (i in 1:R) {
    id <- Rfast2::Sample.int(n, n)
    dvec <- as.vector( crossprod(x[id, ], y) )
    pf <- quadprog::solve.QP( Dmat = Dmat, dvec = dvec, Amat = A, bvec = bvec, meq = px )
    pmse[i] <- pf$value
  }

  ( sum(pmse < mse) + 1 ) / (R + 1)
}
