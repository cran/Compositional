ols.compcomp <- function(y, x, xnew = NULL) {

  py <- dim(y)[2]   ;    px <- dim(x)[2]
  pyx <- py * px    ;    n <- dim(y)[1]

  dvec <-  2 * as.vector( crossprod(x, y) )
  xx <- crossprod(x)
  XX <- matrix(0, pyx, pyx)
  ind <- matrix( 1:pyx, ncol = px, byrow = TRUE)
  for ( i in 1:py )  XX[ ind[i, ], ind[i, ] ] <- xx
  A <- matrix(0, pyx, pyx)
  for ( i in 1:px )  A[i, ind[, i]] <- 1
  A <- t( rbind( A, diag(pyx), -diag(pyx) ) )
  A <- A[, -c( (px + 1): pyx) ]
  bvec <- c( rep(1, px), rep(0, pyx), rep(-1, pyx) )
  f <- quadprog::solve.QP(Dmat = 2 * XX, dvec = dvec, Amat = A, bvec = bvec, meq = px, factorized=FALSE)
  be <- matrix(f$solution, ncol = py)
  mse <- ( sum(y^2) + f$value ) / n

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

  list( mse = mse, be = be, est = est )
}


# ols.compcomp_newer <- function(y, x, xnew = NULL) {

  # py <- dim(y)[2]   ;    px <- dim(x)[2]
  # pyx <- py * px

  # xxinv <- solve( crossprod(x) )
  # m <- numeric( px )
  # dvec <-  - 2 * as.vector( crossprod(x, y) )
  # yy <- sum( diag( crossprod(y) ))

  # ols <- function(be, yy, dvec, m, xxinv) {
    # be <- matrix(be, ncol = py)
    # be <- be / rowSums(be)
    # be1 <- as.vector(be)
    # yy + sum(dvec * be1) + sum( Rfast::mahala(t(be), m, xxinv) )
  # }

  # runtime <- proc.time()
  # mod <- optim( runif(pyx), ols, yy = yy, dvec = dvec, m = m, xxinv = xxinv,
                # method = "L-BFGS-B", lower = rep(0, pyx),
                # upper = rep(1, pyx), control = list(maxit = 10000) )
  # mod <- optim( mod$par, ols, yy = yy, dvec = dvec, m = m, xxinv = xxinv,
                # method = "L-BFGS-B", lower = rep(0, pyx),
                # upper = rep(1, pyx), control = list(maxit = 10000) )
  # runtime <- proc.time() - runtime

  # be <- mod$par
  # be <- matrix(be, ncol = py)
  # be <- be / rowsums(be)


  # if ( is.null( colnames(y) ) ) {
    # colnames(be) <- paste("Y", 1:py, sep = "")
  # } else colnames(be) <- colnames(y)
  # if ( is.null( rownames(y) ) ) {
    # rownames(be) <- paste("X", 1:px, sep = "")
  # } else rownames(be) <- colnames(x)

  # est <- NULL
  # if ( !is.null(xnew) ) {
    # est <- xnew %*% be
  # }

  # list( runtime = runtime, mse = mod$value / dim(y)[1], be = be, est = est )
# }




# ols.compcomp_old <- function(y, x, xnew = NULL) {

  # py <- dim(y)[2]   ;    px <- dim(x)[2]
  # pyx <- py * px

  # ols <- function(be) {
    # be <- matrix(be, ncol = py)
    # be <- be / rowsums(be)
    # mu <- x %*% be
    # sum( (y - mu)^2 )
  # }

  # runtime <- proc.time()
  # mod <- optim( runif(pyx), ols, method = "L-BFGS-B", lower = rep(0, pyx),
                # upper = rep(1, pyx), control = list(maxit = 10000) )
  # mod <- optim( mod$par, ols, method = "L-BFGS-B", lower = rep(0, pyx),
                # upper = rep(1, pyx), control = list(maxit = 10000) )
  # runtime <- proc.time() - runtime

  # be <- mod$par
  # be <- matrix(be, ncol = py)
  # be <- be / rowsums(be)


  # if ( is.null( colnames(y) ) ) {
    # colnames(be) <- paste("Y", 1:py, sep = "")
  # } else colnames(be) <- colnames(y)
  # if ( is.null( rownames(y) ) ) {
    # rownames(be) <- paste("X", 1:px, sep = "")
  # } else rownames(be) <- colnames(x)

  # est <- NULL
  # if ( !is.null(xnew) ) {
    # est <- xnew %*% be
  # }

  # list( runtime = runtime, mse = mod$value / dim(y)[1], be = be, est = est )
# }
