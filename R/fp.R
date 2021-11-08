fp <- function(x, lambda) {

  if ( abs(lambda) > 1e-9 ) {
    y <- x[, -1]^lambda - x[, 1]^lambda
  } else {
    y <- Rfast::Log(x[, -1] / x[, 1])
  }

  y
}
