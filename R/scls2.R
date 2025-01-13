scls2 <- function(y, x, wei = FALSE, xnew = NULL) {

  dx <- length(x)
  px <- numeric(dx)
  xc <- NULL
  for (i in 1:dx)  {
    xc <- cbind(xc, x[[ i ]] / dx )
    px[i] <- dim(x[[ i]] )[2]
  }
  mod <- Compositional::scls(y, xc)
  ini.be <- abs(mod$be)
  px <- c(0, cumsum(px) )

  est <- NULL
  if ( !is.null(xnew) ) {
    est <- 0
    for (i in 1:dx)  est <- est + xnew[[ i ]] %*% ini.be[(px[i] + 1):px[i + 1], ] / dx
  }

  weights <- am <- NULL
  mse <- mod$mse

  if (wei) {
    reg <- function(a, y = y, est = est) {
      b <- c(1, exp(a) )
      b <- b / sum(b)
      b <- round(b, 6)
      ma <- 0
      for ( i in 1:length(b) )  ma <- ma + b[i] * est[[ i ]]
      #- 2 * sum( diag( crossprod(y, ma) ) ) + sum( diag( crossprod(ma) ) )
      - 2 * sum( y * ma ) + sum( ma^2)
    }

    est <- list()
    for (i in 1:dx)  est[[ i ]] <- x[[ i ]] %*% ini.be[(px[i] + 1):px[i + 1], ]
    suppressWarnings( f <- optim( rnorm(dx - 1), reg, y = y, est = est, control = list(maxit = 5000) ) )
    a <- f$par
    b <- c(1, exp(a) )
    b <- b / sum(b)
    am <- round(b, 6)
    mse <- ( sum(y^2) + f$value) / dim(y)[1]
    be <- NULL
    for (i in 1:dx)  be <- rbind( be, am[i] * ini.be[(px[i] + 1):px[i + 1], ] )
    weights <- Rfast::rowsums(be)
    est <- NULL
    if ( !is.null(xnew) ) {
      est <- 0
      for (i in 1:dx) est <- est + xnew[[ i ]] %*% be[(px[i] + 1):px[i + 1], ]
    }
  }

  list(ini.mse = mod$mse, be = ini.be, mse = mse, weights = weights, am = am, est = est)
}
