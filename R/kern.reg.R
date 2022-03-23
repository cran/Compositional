kern.reg <- function(xnew, y, x, h = seq(0.1, 1, length = 10), type = "gauss") {

  y <- as.matrix(y)
  x <- as.matrix(x)
  xnew <- as.matrix(xnew)
  d <- dim(y)[2]
  p <- dim(x)[2]
  n <- dim(y)[1]
  nu <- dim(xnew)[1]
  m <- Rfast::colmeans(x)
  s <- Rfast::colVars(x, std = TRUE)
  x <- t( ( t(x) - m ) / s )  ## standardize the independent variables
  xnew <- t( ( t(xnew) - m ) / s )  ## standardize the x values

  if ( length(h) == 1 ) {
    if (type == "gauss") {
      a1 <- 0.5 * Rfast::dista(xnew, x, square = TRUE)/h^2
    } else  a1 <- Rfast::dista(xnew, x, type = "manhattan" )/h
    z <- exp(-a1)

    ta <- Rfast::rowsums(z)
    mhx <- ( z %*% y) / ta
    mhx[ is.na(mhx) ] <- 0

    if ( is.null(colnames(y)) ) {
      colnames(mhx) <- paste("yhat", 1:d, sep = "" )
    } else  colnames(mhx) <- colnames(y)
    if  ( d == 1 )  mhx <- as.vector(mhx)

  } else {
    len <- length(h)

    if (type == "gauss") {
      a1 <- 0.5 * Rfast::dista(xnew, x, square = TRUE)
      h <- h^2
    } else  a1 <- Rfast::dista(xnew, x, type = "manhattan" )

    if ( d == 1 ) {
      mhx <- matrix(nrow = nu, ncol = len)
      for (i in 1:len) {
        z <- exp( -a1 / h[i] )
        ta <- Rfast::rowsums(z)
        mhx[, i] <- ( z %*% y) / ta
        z <- NULL
      }
      mhx[ is.na(mhx) ] <- 0
      colnames(mhx) <- paste("h=", h, sep = "")

    } else {
      names <- paste("h=", h, sep = "")
      mhx <- sapply(names, function(x) NULL)
      for (i in 1:len) {
        z <- exp( -a1 / h[i] )
        ta <- Rfast::rowsums(z)
        mhx[[ i ]] <- ( z %*% y) / ta
        z <- NULL
      }
    }

  }

  mhx
}

