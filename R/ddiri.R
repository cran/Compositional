################################
#### Dirichlet distribution parameters
#### Tsagris Michail 3/2012
#### mtsagris@yahoo.gr
################################

ddiri <- function(x, a, logged = TRUE) {
  ## x is the compositional data
  ## a is a vector with the parameters

  if ( is.null(nrow(x)) ) {
    x <- x / sum(x)
    f <- lgamma( sum(a) ) - sum( lgamma(a) ) + sum( log(x) * (a - 1) )

  } else {

    x <- as.matrix(x)  ## makes sure x is a matrix
    x <- x / as.vector( Rfast::rowsums(x) )  ## makes sure x is compositional data
    f <- lgamma( sum(a) ) - sum( lgamma(a) ) + as.vector( log(x) %*% (a - 1) )
  }

  if ( logged == TRUE ) {
    f <- f

  } else {
    f <- exp(f)
  }

  f

}

