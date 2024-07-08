tflr <- function(y, x, xnew = NULL) {

  runtime <- proc.time()
  B <- Compositional::scls(y, x)$be
  Dr <- dim(y)[2]   ;    Dp <- dim(x)[2]  ;  Drp <- Dr * Dp
  id <- matrix(1:Drp, ncol = Dr)

  X <- NULL
  for (j in 1:Dr)  X <- cbind(X, x)
  b <- as.vector(B)
  pijk <- Rfast::eachrow(X, b, oper = "*")
  for (j in 1:Dr)  {
    pijk[, id[, j]] <- pijk[, id[, j]] / Rfast::rowsums( pijk[, id[, j]] )
    pijk[ which(is.na(pijk)) ] <- 0
    B[, j] <- Rfast::eachcol.apply( pijk[, id[, j]], y[, j] )
  }
  B <- B / Rfast::rowsums(B)
  yest <- x %*% B
  a <- y * log(yest)
  a[is.infinite(a)] <- NA
  f1 <-  - sum(a, na.rm = TRUE)

  b <- as.vector(B)
  pijk <- Rfast::eachrow(X, b, oper = "*")
  for (j in 1:Dr) {
    pijk[, id[, j]] <- pijk[, id[, j]] / Rfast::rowsums( pijk[, id[, j]] )
    pijk[ which(is.na(pijk)) ] <- 0
    B[, j] <- Rfast::eachcol.apply( pijk[, id[, j]], y[, j] )
  }
  B <- B / Rfast::rowsums(B)
  yest <- x %*% B
  a <- y * log(yest)
  a[is.infinite(a)] <- NA
  f2 <-  - sum(a, na.rm = TRUE)
  i <- 2

  while(f1 - f2 > 1e-6) {
    i <- i + 1
    f1 <- f2
    b <- as.vector(B)
    pijk <- Rfast::eachrow(X, b, oper = "*")
    for (j in 1:Dr) {
      pijk[, id[, j]] <- pijk[, id[, j]] / Rfast::rowsums( pijk[, id[, j]] )
      pijk[ which(is.na(pijk)) ] <- 0
      B[, j] <- Rfast::eachcol.apply( pijk[, id[, j]], y[, j] )
    }
    B <- B / Rfast::rowsums(B)
    yest <- x %*% B
    a <- y * log(yest)
    a[is.infinite(a)] <- NA
    f2 <-  - sum(a, na.rm = TRUE)
  }
  runtime <- proc.time() - runtime

  if ( is.null( colnames(y) ) )  {
    colnames(B) <- paste("Y", 1:Dr, sep = "")
  } else colnames(B) <- colnames(y)

  if ( is.null( colnames(x) ) )  {
    rownames(B) <- paste("X", 1:Dp, sep = "")
  } else rownames(B) <- colnames(x)

  kl <- sum( y * log(y), na.rm = TRUE ) + f2
  est <- NULL
  if ( !is.null(xnew) )  est <- xnew %*% B

  list(runtime = runtime, iters = i, kl = kl, be = B, est = est)
}








# .tflr <- function(y, x, xnew = NULL) {
#   py <- dim(y)[2]   ;    px <- dim(x)[2]
#   pyx <- py * px    ;    n <- dim(y)[1]
#   n <- dim(y)[1]    ;    npy <- n * py
#
#   X <- matrix(0, npy, pyx)
#   indr <- matrix( 1:npy, ncol = py )
#   indc <- matrix( 1:pyx, ncol = py )
#   for ( i in 1:py )  X[ indr[, i], indc[, i] ] <- x
#   Y <- as.vector(y)
#
#   ind <- matrix( 1:pyx, ncol = px, byrow = TRUE )
#   A <- matrix(0, pyx, pyx)
#   for ( i in 1:px )  A[i, ind[, i]] <- 1
#   A <- t( rbind( A, diag(pyx), -diag(pyx) ) )
#   A <- A[, -c( (px + 1): pyx) ]
#   bvec <- c( rep(1, px), rep(0, pyx), rep(-1, pyx) )
#   A <- t(A)
#
#   a <- goric::orglm(Y ~ X - 1, family = quasibinomial(link="identity"),
#                     data = data.frame(Y = Y, X = X), constr = A, rhs = bvec, nec = px)
#   be <- matrix( abs(coef(a)), ncol = py)
#   est <- x %*% be
#   kl <- sum(y * log(y / est), na.rm = TRUE)
#
#   if ( is.null( colnames(y) ) ) {
#     colnames(be) <- paste("Y", 1:py, sep = "")
#   } else colnames(be) <- colnames(y)
#   if ( is.null( colnames(x) ) ) {
#     rownames(be) <- paste("X", 1:px, sep = "")
#   } else rownames(be) <- colnames(x)
#
#   est <- NULL
#   if ( !is.null(xnew) ) {
#     est <- xnew %*% be
#   }
#
#   list( kl = kl, be = be, est = est )
# }






# tflr <- function(y, x, xnew = NULL) {
#
#   runtime <- proc.time()
#   be <- codalm::codalm(y, x)
#   runtime <- proc.time() - runtime
#   d <- dim(y)[2]
#   p <- dim(x)[2]
#
#   if ( is.null( colnames(y) ) )  {
#     colnames(be) <- paste("Y", 1:d, sep = "")
#   } else colnames(be) <- colnames(y)
#
#   if ( is.null( rownames(y) ) )  {
#     rownames(be) <- paste("X", 1:p, sep = "")
#   } else rownames(be) <- colnames(x)
#
#   yhat <- x %*% be
#   kl <- sum( y * log(y / yhat), na.rm = TRUE )
#   est <- NULL
#   if ( !is.null(xnew) )  est <- xnew %*% be
#
#   list(runtime = runtime, kl = kl, be = be, est = est)
# }
