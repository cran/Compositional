ols.compcomp <- function(y, x, rs = 5, xnew = NULL) {

  con <- function(be){
    be <- matrix(be, byrow = TRUE, ncol = py)
    f <- rowSums(be) - 1
    list(ceq = f, c = NULL)
  }

  ols <- function(be) {
    be <- matrix(be, byrow = TRUE, ncol = py)
    mu <- x %*% be
    sum( (y - mu)^2 )
  }

  py <- dim(y)[2]   ;    px <- dim(x)[2]
  pyx <- py * px

  sse <- numeric(rs)
  bers <- matrix(nrow = pyx, ncol = rs)

  runtime <- proc.time()
  for (i in 1:rs) {
    f1 <- NlcOptim::solnl( X = runif(pyx), ols, con, lb = rep(0, pyx), ub = rep(1, pyx) )
    f2 <- NlcOptim::solnl( f1$par, ols, con, lb = rep(0, pyx), ub = rep(1, pyx) )
    while ( f1$fn - f2$fn > 1e-04 ) {
      f1 <- f2
      f1 <- NlcOptim::solnl( f2$par, ols, con, lb = rep(0, pyx), ub = rep(1, pyx) )
      f2 <- NlcOptim::solnl( f1$par, ols, con, lb = rep(0, pyx), ub = rep(1, pyx) )
    }
    sse[i] <- f2$fn
    bers[, i] <- f2$par
  }
  runtime <- proc.time() - runtime

  k <- which.min(sse)
  be <- matrix(bers[, k], byrow = TRUE, ncol = py)
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

  list( runtime = runtime, mse = sse[k]/dim(y)[1], be = be, est = est )
}
