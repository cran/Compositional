mlr <- function(x) {
  d <- dim(x)[2] - 1
  x <- t(x)
  t( log( x[1:d, ] / ( 1 - Rfast::colCumSums(x)[1:d, ] ) ) )
}
