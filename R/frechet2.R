frechet2 <- function(x, di, a, k) {
  n <- dim(di)[1]   ;   p <- length(k)
  d <- dim(x)[2]
  knam <- paste("k=", k, sep = "")
  m <- sapply(knam, function(x) NULL)
  denom <- 1:max(k)
  m1 <- matrix(nrow = n, ncol = d * length(k) )
  di <- t(di)
  
  if ( abs(a) < 1e-9 ) {
  
    for ( i in 1:n ) {
      lx <- Rfast::colCumSums( Rfast::Log( x[ di[, i], ] ) ) / denom
      esk <- exp( lx )
      est <- esk/Rfast::rowsums(esk)
      m1[i, ] <- as.vector( t(est[-1, ]) )
    }

  } else {
    inva <- 1/a
    for ( i in 1:n ) {
      xa <- x[ di[, i], ]^a
      z <- xa / Rfast::rowsums(xa)
      esk <- Rfast::colCumSums(z)^inva / denom
      est <- esk/Rfast::rowsums(esk)
      m1[i, ] <- as.vector( t( est[k, , drop = FALSE] ) )
    }

  }

  ind <- matrix( 1:dim(m1)[2], ncol = p )
  for ( j in 1:p)  m[[ j ]] <- m1[, ind[, j] ]
  m
}

