ulc.reg <- function(y, x, z = NULL, xnew = NULL, znew = NULL) {

  if ( is.null(z) ) {
    x <- log(x)
    x <- cbind(1, x)
    dm <- dim(x)
    n <- dm[1]  ;  p <- dm[2]
    xxs <- solve( crossprod(x) )
    be <- as.vector( xxs %*% crossprod(x, y) )
    e <- y - x %*% be
    va <- sum(e^2) / (n - p)
    covbe <- xxs * va
    nama <- colnames(x)
    if ( is.null(nama) )  nama <- c( "constant", paste("X", 1:(p-1), sep = "") )
    if ( nama[1] == "" )  nama[1] <- "constant"
    names(be) <- nama
    colnames(covbe) <- rownames(covbe) <- nama
    est <- NULL
    if ( !is.null(xnew) )  est <- cbind(1, log(xnew) ) %*% be

  } else {
    z <- model.matrix(y~., data = as.data.frame(z) )[, -1]
    X <- cbind(1, log(x), z )
    dm <- dim(X)
    n <- dm[1]  ;  p <- dm[2]
    xxs <- solve( crossprod(X) )
    be <- as.vector( xxs %*% crossprod(X, y) )
    e <- y - X %*% be
    va <- sum(e^2) / (n - p)
    covbe <- xxs * va
    nama <- colnames(X)
    d1 <- dim(x)[2]
    d2 <- p - d1 - 1
    if ( is.null(nama) )  nama <- c( "constant", paste("X", 1:d1, sep = ""),
                                     paste("Z", 1:d2, sep = "") )
    if ( nama[1] == "" )  nama[1] <- "constant"
    names(be) <- nama
    colnames(covbe) <- rownames(covbe) <- nama
    est <- NULL
    if ( !is.null(xnew)  &  !is.null(znew) ) {
      znew <- model.matrix(~., data.frame(znew))[, -1]
      est <- cbind(1, log(xnew), znew) %*% be
    }
  }

  list(be = be, covbe = covbe, va = va, residuals = e, est = est)
}







