################################
#### Dirichlet distribution parameters
#### Tsagris Michail 3/2012
#### mtsagris@yahoo.gr
################################

ddiri <- function(x, a, logged = TRUE) {
  ## x is the compositional data
  ## a is a vector with the parameters
  if ( is.null(dim(x)[1]) ) {
    f <- lgamma( sum(a) ) - sum( lgamma(a) ) + sum( log(x) * (a - 1) )
  } else  f <- lgamma( sum(a) ) - sum( lgamma(a) ) + as.vector( log(x) %*% (a - 1) )

  if ( logged ) {
    f <- f
  } else   f <- exp(f)

  f
}

