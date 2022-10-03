################################
#### Dirichlet distribution parameters
#### Tsagris Michail 3/2012
#### mtsagris@yahoo.gr
################################
dgendiri <- function(x, a, b, logged = TRUE) {
  ## x is the compositional data
  ## a is a vector with the parameters
  if ( is.null(dim(x)[1]) ) {
    f <- lgamma( sum(a) ) - sum( lgamma(a) ) + sum( a * log(b) ) + sum( (a - 1) * log(x) ) -
      sum(a) * log( sum(b * x) )
  } else  f <- lgamma( sum(a) ) - sum( lgamma(a) ) + sum( a * log(b) ) +
      as.vector( log(x) %*% (a - 1) ) - sum(a) * as.vector( log(x %*% b) )

  if ( logged ) {
    f <- f
  } else   f <- exp(f)

  f
}

