scrq <- function(y, x, xnew = NULL) {

  py <- dim(y)[2]   ;    px <- dim(x)[2]
  pyx <- py * px    ;    n <- dim(y)[1]
  n <- dim(y)[1]    ;    npy <- n * py

  X <- matrix(0, npy, pyx)
  indr <- matrix( 1:npy, ncol = py )
  indc <- matrix( 1:pyx, ncol = py )
  for ( i in 1:py )  X[ indr[, i], indc[, i] ] <- x
  Y <- as.vector(y)

  R <- NULL
  for (i in 1:py)  R <- cbind(R, diag(px))
  R <- rbind( R, -R, diag(pyx), -diag(pyx) )
  r <- c( rep(1, px), rep(-1, px), rep(0, pyx), rep(-1, pyx) )

  a <- quantreg::rq(Y ~ X - 1, data = data.frame(Y = Y, X = X), method = "fnc", R = R, r = r)
  be <- matrix(coef(a), ncol = py)
  mlad <- sum( abs (y - x %*% be) ) / n

  if ( is.null( colnames(y) ) ) {
    colnames(be) <- paste("Y", 1:py, sep = "")
  } else colnames(be) <- colnames(y)
  if ( is.null( rownames(y) ) ) {
    rownames(be) <- paste("X", 1:px, sep = "")
  } else rownames(be) <- colnames(x)

  est <- NULL
  if ( !is.null(xnew) ) {
    est <- xnew %*% be
  }

  list( mlad = mlad, be = be, est = est )
}
