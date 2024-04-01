tflr <- function(y, x, xnew = NULL) {

  runtime <- proc.time()
  B <- Compositional::scls(y, x)$be
  Dr <- dim(y)[2]   ;    Dp <- dim(x)[2]  ;  Drp <- Dr * Dp
  id <- matrix(1:Drp, ncol = Dr)

  X <- NULL
  for (j in 1:Dr)  X <- cbind(X, x)
  b <- as.vector(B)
  pijk <- Rfast::eachrow(X, b, oper = "*")
  for (j in 1:Dr)  pijk[, id[, j]] <- pijk[, id[, j]] / Rfast::rowsums( pijk[, id[, j]] )
  pijk[ which(is.na(pijk)) ] <- 0
  for (j in 1:Dr)  B[, j] <- Rfast::eachcol.apply( pijk[, id[, j]], y[, j] )
  B <- B / Rfast::rowsums(B)
  yest <- x %*% B
  a <- y * log(yest)
  a[is.infinite(a)] <- NA
  f1 <-  - sum(a, na.rm = TRUE)

  b <- as.vector(B)
  pijk <- Rfast::eachrow(X, b, oper = "*")
  for (j in 1:Dr)  pijk[, id[, j]] <- pijk[, id[, j]] / Rfast::rowsums( pijk[, id[, j]] )
  pijk[ which(is.na(pijk)) ] <- 0
  for (j in 1:Dr)  B[, j] <- Rfast::eachcol.apply( pijk[, id[, j]], y[, j] )
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
    for (j in 1:Dr)  pijk[, id[, j]] <- pijk[, id[, j]] / Rfast::rowsums( pijk[, id[, j]] )
    pijk[ which(is.na(pijk)) ] <- 0
    for (j in 1:Dr)  B[, j] <- Rfast::eachcol.apply( pijk[, id[, j]], y[, j] )
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

  if ( is.null( rownames(y) ) )  {
    rownames(B) <- paste("X", 1:Dp, sep = "")
  } else rownames(B) <- colnames(x)

  kl <- sum( y * log(y), na.rm = TRUE ) + f2
  est <- NULL
  if ( !is.null(xnew) )  est <- xnew %*% B

  list(runtime = runtime, iters = i, kl = kl, be = B, est = est)
}



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
