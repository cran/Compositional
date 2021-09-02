pivotinv <- function(y) {
  x <- cbind(0, y)
  p <- dim(x)[2]
  y <- t(y)
  x[, 1] <- sqrt( (p - 1 ) / p ) * y[1, ]
  for ( j in 2:c(p - 1) ) {
    k <- 1:c(j - 1)
    x[, j] <-  - Rfast::colsums( 1 / sqrt( (p - k + 1) * (p - k) ) * y[k, , drop = FALSE] ) + 
              sqrt( (p - j) / (p - j + 1) ) * y[j, ]
  }
  k <- 1:c(p - 1)
  x[, p] <-  - Rfast::colsums( 1 / sqrt( (p - k + 1) * (p - k) ) * y[k, , drop = FALSE] ) 
  x <- exp(x)
  x / Rfast::rowsums(x)
}