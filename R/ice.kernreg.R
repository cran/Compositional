ice.kernreg <- function(y, x, h, type = "gauss", k = 1, frac = 0.1) {

  x <- as.matrix(x)
  n <- dim(x)[1]
  nu <- ceiling( frac * n )
  xsel <- sort( Rfast2::Sample(x[, k], nu) )
  est <- matrix(NA, n, nu)

  for ( i in 1:nu ) {
    X <- x
    X[, k] <- xsel[i]
    est[, i] <- Compositional::kern.reg(X, y, x, h = h, type = "gauss")
  }

  nam <- colnames(x)[k]
  if ( is.null(nam) )  nam <- paste("Variable ", k, sep = "")

  est <- est[, -1] - est[, 1]
  plot( xsel[-1], est[1, ], type = "l", xlab = nam, ylab = "Centered fitted values",
        cex.lab = 1.3, cex.axis = 1.3, ylim = c( min(est), max(est) ) )
  abline(v = seq( min(xsel), max(xsel), length = 10 ), col = "lightgrey", lty = 2)
  abline(h = seq(min(est), max(est), length = 10), col = "lightgrey", lty = 2)
  for (i in 2:n )  lines(xsel[-1], est[i, ])
  m <- Rfast::colmeans(est)
  lines(xsel[-1], m, col = 4, lwd = 3)
}
