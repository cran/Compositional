alfa <- function(x, a, h = TRUE) {
  ## x contains the compositional data
  ## a is the power parameter, usually between -1 and 1
  ## if h is TRUE the multiplication with the Helmert matrix takes place
  x <- as.matrix(x)
  D <- dim(x)[2] ## number of components
  if ( D == 1 )   x <- t(x)
  D <- dim(x)[2] ## number of components
  
  if ( abs(a) > 1e-9 ) {
    z <- x^a
    ta <- Rfast::rowsums(z)
    z <- D / a * z / ta - 1/a
    sa <- sum( log(ta) )
  } else {  ## if a = 0 the ilr is calculated
    xa <- Rfast::Log(x)
    z <- xa - Rfast::rowmeans( xa )   ## this is the clr
    sa <- dim(x)[1] * log(D)
  }

  if ( h ) {
    aff <- tcrossprod(z, helm( D ) ) ## multiply by the Helmert sub-matrix
    res <- list(sa = sa, aff = aff)
  } else  res <- list(sa = sa, aff = z)

  res
}
