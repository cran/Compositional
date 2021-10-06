ice.pprcomp <- function(y, x, nterms = 3, k = 1, type = "alr") {
  
 if ( type == "alr" ) {
    x <- Compositional::alr(x)
  } else  x <- Rfast::Log(x)

  x <- as.data.frame(x)
  nam <- colnames(x)
  n <- dim(x)[1]   ;   p <- dim(x)[2]
  if ( is.null(nam) )  colnames(x) <- paste("X", 1:p, sep = "")
  nam <- colnames(x)
  xsel <- sort( x[, k] )  
  X <- x
  est <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      X[i, k] <- xsel[j]    
      mod <- ppr(y ~., data = X, nterms = nterms)
      est[i, j] <- predict(mod, newdata = X)[1]
    }
  }
  est <- est[, -1] - est[, 1]
  plot( xsel, est[, 1], type = "l", xlab = nam[k], ylab = "Centered fitted values",
        cex.lab = 1.3, cex.axis = 1.3, ylim = c( min(est), max(est) ) )
  for (j in 2:c(n - 1) )  lines(xsel, est[, j])
  m <- Rfast::rowmeans(est)
  lines(xsel, m, col = 4, lwd = 3)
  a <- loess( m ~ xsel )
  lines(xsel, fitted(a), col = 5, lwd = 2)
}
  