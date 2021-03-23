zeroreplace <- function(x, a = 0.65, delta = NULL, type = "multiplicative") {

  d <- dim(x)[2]
  if ( is.null(delta) ) {
    delta <- a * min(x[x > 0] )
  }
  pou <- which(x == 0, arr.ind = TRUE )
  len <- unique(pou[, 1])

  if ( type == "multiplicative" ) {
    for (i in len) {
      co <- which( x[i, ] == 0)
      z <- length(co)
      x[i, co] <- delta
      x[i, -co] <- ( 1 - z * delta ) * x[i, -co]
    }

  } else if ( type == "additive" ) {
    for (i in len) {
      co <- which( x[i, ] == 0)
      z <- length(co)
      x[i, co] <- delta * (z + 1) * (d - z) / d^2
      x[i, -co] <- x[i, -co] - delta * (z + 1) * z / d^2
    }

  } else if ( type == "simple" ) {
    for (i in len) {
      co <- which( x[i, ] == 0)
      z <- length(co)
      denom <- 1 + z * delta
      x[i, co] <- delta / denom
      x[i, -co] <- x[i, -co] / denom
    }
    x[len, ] <- x[len, ]/Rfast::rowsums(x[len, ])

  }

  x
}
