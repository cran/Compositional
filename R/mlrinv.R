mlrinv <- function(y) {
  d <- dim(y)[2]
  x <- t( exp(y) )
  a <- Rfast::colCumProds(1 + x)
  cbind( t( x / a), 1 / a[d, ] )
}
