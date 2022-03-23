ice.pprcomp <- function(model, x, k = 1, frac = 0.1, type = "log") {

  nam <- colnames(x)
  if ( type == "alr" ) {
    x <- Compositional::alr(x)
    colnames(x) <- nam[-1]
  } else {
    x <- Rfast::Log(x)
    colnames(x) <- nam
  }

  mod <- model$mod
  dm <- dim(x)
  n <- dm[1]  ;  p <- dm[2]
  x <- as.data.frame(x)
  nu <- ceiling( frac * n )
  xsel <- sort( sample(x[, k], nu) )
  est <- matrix(NA, n, nu)

  for ( i in 1:nu ) {
    X <- x
    X[, k] <- xsel[i]
    est[, i] <- predict(mod, newdata = X)
  }

  nam <- colnames(x)[k]
  if ( is.null(nam) )  nam <- paste("Variable ", k, sep = "")

  est <- est[-1, ] - est[1, ]
  plot( xsel, est[1, ], type = "l", xlab = nam, ylab = "Centered fitted values",
        cex.lab = 1.3, cex.axis = 1.3, ylim = c( min(est), max(est) ) )
  abline(v = seq( min(xsel), max(xsel), length = 10 ), col = "lightgrey", lty = 2)
  abline(h = seq(min(est), max(est), length = 10), col = "lightgrey", lty = 2)
  for (i in 2:c(n - 1) )  lines(xsel, est[i, ])
  m <- Rfast::colmeans(est)
  lines(xsel, m, col = 4, lwd = 3)
}
