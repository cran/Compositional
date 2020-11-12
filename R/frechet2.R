################################
#### Frechet mean
#### Tsagris Michail 5/2013
#### References: Tsagris, M. T., Preston, S., and Wood, A. T. A. (2011).
#### A data-based power transformation for
#### compositional data. In Proceedings of the 4rth Compositional Data Analysis Workshop, Girona, Spain.
#### mtsagris@yahoo.gr
################################
frechet2 <- function(x, di, a, k1) {
  dm <- dim(di)
  n <- dm[1]   ;   p <- dm[2]
  d <- dim(x)[2]
  knam <- paste("k=", p - 1, sep = "")
  m <- sapply(knam, function(x) NULL)
  denom <- 1:p
  m1 <- matrix(nrow = n, ncol = d * (p - 1) )
  apo <- 1:k1
  if ( a == 0 ) {
    for ( i in 1:n ) {
      lx <- Rfast::colCumSums( Rfast::Log( x[ di[i, ], ] ) ) / denom
      esk <- exp( lx[-apo, , drop = FALSE] )
      est <- esk/Rfast::rowsums(esk)
      m1[i, ] <- as.vector( t(est) )
    }

  } else {
    inva <- 1/a
    for ( i in 1:n ) {
      xa <- x[ di[i, ], ]^a
      z <- xa / Rfast::rowsums(xa)
      esk <- Rfast::colCumSums(z)^inva / denom
      est <- esk/Rfast::rowsums(esk)
      m1[i, ] <- as.vector( t( est[-apo, , drop = FALSE] ) )
    }

  }

  ind <- matrix( 1:dim(m1)[2], ncol = p - 1)
  for ( j in 1:(p - 1) )  m[[ j ]] <- m1[, ind[, j] ]
  m
}
