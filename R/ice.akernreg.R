ice.akernreg <- function(y, x, a, h, type = "gauss", ind = 1, frac = 0.1, qpos = 0.9) {

  dm <- dim(y)
  n <- dm[1]  ;  d <- dm[2]
  nu <- ceiling(frac * n)
  x <- as.matrix(x)
  xsel <- sort( Rfast2::Sample(x[, ind], nu) )
  est <- matrix(NA, nu, d)
  for (i in 1:nu) {
    X <- x
    X[, ind] <- xsel[i]
    est[i, ] <- Rfast::colmeans( Compositional::akern.reg(X, y, x, a, h, type)[[ 1 ]][[ 1 ]] )
  }

  namx <- colnames(x)[ind]
  if ( is.null(namx) )  namx <- paste("Variable ", ind, sep = "")

  ##png(filename = "ice.png", width = 5000, height = 4000, res = 600)

  par(mar = c(8, 5.5, 8, 9.5), xpd = TRUE)
  plot( xsel, est[, 1], type = "l", xlab = namx, ylab = "Fitted proportions",
        xaxt = "n", yaxt = "n", cex.lab = 1.2, cex.axis = 1.2,
        ylim = c( min(est), max(est) ), lwd = 2, bty = "n" )

  v <- seq( min(xsel), max(xsel), length = 10 )
  h <- seq( min(est), max(est), length = 10 )
  mtext( text = round(h, 2), side = 2, at = h, las = 2, font = 2, line = 0.2 )
  mtext( text = round(v, 2), side = 1, at = v, las = 2, font = 2, line = 0.2 )

  for (i in 1:10) {
    lines(rep(v[i], 10), h, col = "lightgrey", lty = 2)
    lines(v, rep(h[i], 10), col = "lightgrey", lty = 2)
  }
  for ( i in 1:d ) lines(xsel, est[, i], col = i, lwd = 2)

  namy <- colnames(y)
  if ( is.null(namy) )  namy <- paste("Comp. ", 1:d, sep = "")

  legend( x = quantile(x[, ind], qpos), y = max(est), legend = namy,
          xpd = TRUE, bty = "n", title = "Components", col = 1:d,
          lwd = rep(2, d), lty = rep(1, d), text.col = 1:d )

  ##dev.off()

}
