zeroreplace <- function(x, a = 2/3) {
  d <- dim(x)[2]
  delta <- a * min(x[x > 0] )
  pou <- which(x == 0, arr.ind = TRUE )
  len <- unique(pou[, 1])
  for (i in len) {
    co <- which( x[i, ] == 0)
    z <- length(co)
    x[i, co] <- delta * (z + 1) * (d - z) / d^2
    x[i, -co] <- x[i, -co] - delta * (z + 1) * z / d^2
  }
  x
}
