tflr2 <- function(y, x, wei = FALSE, xnew = NULL) {

  dx <- length(x)
  px <- numeric(dx)
  xc <- NULL
  for (i in 1:dx)  {
    xc <- cbind(xc, x[[ i ]] / dx )
    px[i] <- dim(x[[ i]] )[2]
  }
  mod <- Compositional::tflr(y, xc)
  ini.be <- mod$be
  px <- c(0, cumsum(px) )

  est <- NULL
  if ( !is.null(xnew) ) {
    est <- 0
    for (i in 1:dx)  est <- est + xnew[[ i ]] %*% ini.be[(px[i] + 1):px[i + 1], ] / dx
  }

  weights <- am <- NULL
  kl <- mod$kl

  if (wei) {
    reg <- function(a, y = y, est = est) {
      b <- c(1, exp(a) )
      b <- b / sum(b)
      b <- round(b, 6)
      ma <- 0
      for ( i in 1:length(b) )  ma <- ma + b[i] * est[[ i ]]
      f <- y * log(est)
      f[is.infinite(f)] <- NA
      sum(f)
    }

    est <- list()
    for (i in 1:dx)  est[[ i ]] <- x[[ i ]] %*% ini.be[(px[i] + 1):px[i + 1], ]
    suppressWarnings( f <- optim( rnorm(dx - 1), reg, y = y, est = est, control = list(maxit = 5000) ) )
    a <- f$par
    b <- c(1, exp(a) )
    b <- b / sum(b)
    am <- round(b, 6)
    kl <- sum(y * log(y), na.rm = FALSE) - f$value
    be <- NULL
    for (i in 1:dx)  be <- rbind(be, am[i] * ini.be[(px[i] + 1):px[i + 1], ] )
    weights <- Rfast::rowsums(be)
    est <- NULL
    if ( !is.null(xnew) ) {
      est <- 0
      for (i in 1:dx) est <- est + xnew[[ i ]] %*% be[(px[i] + 1):px[i + 1], ]
    }
  }

  list(ini.kl = mod$kl, be = ini.be, kl = kl, weights = weights, am = am, est = est)
}

