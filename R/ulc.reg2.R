ulc.reg2 <- function(y, x, z = NULL, xnew = NULL, znew = NULL) {

  if ( is.null(z) ) {
    X <- NULL
    p1 <- 0
    for ( i in 1:length(x) ) {
      X <- cbind(X, x[[ i ]] )
      p1 <- c( p1, dim(X)[2] )
    }
    x <- cbind(1, log(X) )
    dm <- dim(x)
    n <- dm[1]  ;  p <- dm[2]

    xxs <- solve( crossprod(x) )
    be <- xxs %*% crossprod(x, y)
    e <- y - x %*% be
    va <- sum(e^2) / (n - p)
    covbe <- xxs * va
    nama <- colnames(x)
    if ( is.null(nama) ) {
      nama <- "constant"
      p1 <- diff(p1)
      for ( i in 1:3 )  nama <- c(nama, paste("X", i, ".", 1:p1[i], sep = "") )
    }
    names(be) <- nama
    colnames(covbe) <- rownames(covbe) <- nama
    est <- NULL
    if ( !is.null(xnew) )  est <- cbind(1, log(xnew) ) %*% be

  } else {

    z <- model.matrix(y~., data = as.data.frame(z) )[, -1, drop = FALSE]
    d2 <- dim(z)[2]
    X <- NULL
    p1 <- 0
    for ( i in 1:length(x) ) {
      X <- cbind(X, x[[ i ]] )
      p1 <- c( p1, dim(X)[2] )
    }
    x <- cbind( 1, z, log(X) )
    dm <- dim(x)
    n <- dm[1]  ;  p <- dm[2]

    xxs <- solve( crossprod(x) )
    be <- xxs %*% crossprod(x, y)
    e <- y - x %*% be
    va <- sum(e^2) / (n - p)
    covbe <- xxs * va
    nama <- colnames(x)

    if ( is.null(nama) ) {
      nama <- "constant"
      p1 <- diff(p1)
      for ( i in 1:length(p1) )  nama <- c(nama, paste("X", i, ".", 1:p1[i], sep = "") )
      nama <- c(paste("Z", 1:d2, sep = ""), nama )
    }
    est <- NULL
    if ( !is.null(xnew)  &  !is.null(znew) ) {
      znew <- model.matrix( ~., data.frame(znew) )
      Xnew <- NULL
      for ( i in 1:length(xnew) )  Xnew <- cbind(Xnew, xnew[[ i ]] )
      est <- cbind(znew, log(Xnew) ) %*% be
    }

  }

  list(be = be, covbe = covbe, va = va, residuals = e, est = est)
}
