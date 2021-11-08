choose.pc <- function(x, plot = TRUE) {
  ## x contains the data
  n <- dim(x)[1]
  runtime <- proc.time()
  #x <- Rfast::standardise(x, center = TRUE, scale = FALSE)  ## center the matrix
  A <- svd(x)
  u <- A$u
  d <- A$d
  v <- t( A$v )
  p <- length(d)
  press <- numeric(p)

  y <- 0
  for (i in 1:p) {
    y <- y + u[, i, drop = FALSE] %*% ( d[i] * v[i, , drop = FALSE] )
    press[i] <- sum( (y - x)^2 )  ## calculation of the PRESS
  }
  runtime <- proc.time() - runtime
  
  if ( plot ) {
    plot(press, type = "b", pch = 9, xlab = "Number of components",
         ylab = "Reconstruction error", lwd = 2, cex.lab = 1.2, cex.axis = 1.2)
  }
  
  val <- d^2 / (n - 1)
  cumprop <- cumsum(val) / sum(val)
  diffa <- diff( c(cumprop, 1) )
  list(values = val, cumprop = cumprop, differences = diffa, press = sqrt(press),
       runtime = runtime)
}
