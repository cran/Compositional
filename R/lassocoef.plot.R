lassocoef.plot <- function(lasso, lambda = TRUE) {
  mod <- lasso$mod
  p <- dim( mod$beta[[ 1 ]] )[1]
  d <- length(mod$beta)
  be <- matrix( nrow = p, ncol = length(mod$lambda) )
  for (i in 1:p) {
    be[i, ] <- mod$beta[[ 1 ]][i, ]^2
    for (j in 2:d) {
      be[i, ] <- be[i, ] + mod$beta[[ j ]][i, ]^2
    }
  }
  be <- sqrt(be)

  if ( lambda ) {
    d <- dim(be)[2]
    be <- be[, d:1]
    Lone <- log( mod$lambda )
    plot(Lone, be[1, ], xlab = expression( paste(Log(lambda) ) ),
         ylab = expression( paste("Coefficients  ", L[2], "-norm") ),
         cex.lab = 1.2, cex.axis = 1.2, type = "l", lwd = 2, ylim = c( min(be), max(be) ) )
    abline(v = seq(0, max(Lone), length = 10), col = "lightgrey", lty = 2)
    abline(h = seq(0, max(be), length = 10), col = "lightgrey", lty = 2)
    for ( i in 1:p )  lines( Lone, be[i, ], col = i, lwd = 2)
  } else {
    Lone <- numeric( dim(be)[2] )
    for (j in 1:d)  Lone <- Lone + Rfast::colsums( abs( as.matrix( mod$beta[[ j ]] ) ) )
    plot(Lone, be[1, ], xlab = expression( paste(L[1], "-norm") ),
         ylab = expression( paste("Coefficients  ", L[2], "-norm") ),
         cex.lab = 1.2, cex.axis = 1.2, type = "l", lwd = 2, ylim = c( min(be), max(be) ) )
    abline(v = seq(0, max(Lone), length = 10), col = "lightgrey", lty = 2)
    abline(h = seq(0, max(be), length = 10), col = "lightgrey", lty = 2)
    for ( i in 1:p )  lines( Lone, be[i, ], col = i, lwd = 2)
  }
}
