bc <- function(x, lambda) {
  if ( lambda != 0 ) {
    y <- ( ( x[, -1] / x[, 1] )^lambda - 1 ) / lambda
  } else y <- Log( x[, -1]/x[, 1] )
  y
}
