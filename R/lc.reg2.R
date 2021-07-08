lc.reg2 <- function(y, x, z = NULL, xnew = NULL, znew = NULL) {

  if ( is.null(z) ) {
    X <- NULL
    p1 <- 0
    for ( i in 1:length(x) ) {
      X <- cbind(X, x[[ i ]] )
      p1 <- c( p1, dim(X)[2] )
    }
    R <- matrix( 0, dim(X)[2], length(p1) - 1 )
    for ( j in 1: (length(p1) - 1) )  R[(p1[j] + 1) : p1[j + 1], j] <- 1
    R <- rbind(0, R)
    x <- cbind(1, log(X) )
    dm <- dim(x)
    n <- dm[1]  ;  p <- dm[2]
    ca <- numeric(p)

    xxs <- solve( crossprod(x) )
    bols <- xxs %*% crossprod(x, y)
    com <- xxs %*% R %*% solve( t(R) %*% xxs %*% R )
    be <- as.vector( bols - com %*% t(R) %*% bols ) - ca
    e <- y - x %*% be
    va <- sum(e^2) / (n - p)
    covbe <- ( xxs - com %*% t(R) %*% xxs ) * va
    nama <- colnames(x)
    if ( is.null(nama) ) {
      nama <- "constant"
      p1 <- diff(p1)
      for ( i in 1:length(p1) )  nama <- c(nama, paste("X", i, ".", 1:p1[i], sep = "") )
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
    R <- matrix( 0, dim(X)[2], length(p1) - 1 )
    for ( j in 1: (length(p1) - 1) )  R[(p1[j] + 1) : p1[j + 1], j] <- 1
    R <- rbind(0, R)
    R <- rbind( matrix(0, d2, dim(R)[2]), R )
    x <- cbind( 1, z, log(X) )
    dm <- dim(x)
    n <- dm[1]  ;  p <- dm[2]
    ca <- numeric(p)

    xxs <- solve( crossprod(x) )
    bols <- xxs %*% crossprod(x, y)
    com <- xxs %*% R %*% solve( t(R) %*% xxs %*% R )
    be <- as.vector( bols - com %*% t(R) %*% bols ) - ca
    e <- y - x %*% be
    va <- sum(e^2) / (n - p)
    covbe <- ( xxs - com %*% t(R) %*% xxs ) * va
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
