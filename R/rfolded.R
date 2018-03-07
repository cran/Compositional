rfolded <- function(n, mu, su, a) {

  z <- Rfast::rmvnorm(n, mu, su)
  D <- dim(z)[2]
  z1 <- z %*% helm(D + 1)
  mi <- abs(a) * Rfast::rowMins(z1, value = TRUE)
  ina <- which(mi <  - 1)
  if ( length(ina) > 0 ) {
    p <- 1 / abs( mi[ina] )
    w <- p^2 * z1[ina, ]
    x1 <- rbind( z1[-ina, ], w )
  } else  x1 <- z1

  if ( abs(a) < 1e-10 )  {  ## if alpha is almost zero make it zero
    z <- exp(x1)
  } else  z <- (a * x1 + 1) ^ (1/a)
    z <- z / Rfast::rowsums(z)
  z
}
