bc <- function(x, lambda) {

  if ( abs(lambda) > 1e-9 ) {
    y <- ( ( x[, -1] / x[, 1] )^lambda - 1 ) / lambda
  } else y <- Log( x[, -1]/x[, 1] )

  y
}
