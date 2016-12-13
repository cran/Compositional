################################
#### Frechet mean
#### Tsagris Michail 5/2013
#### References: Tsagris, M. T., Preston, S., and Wood, A. T. A. (2011).
#### A data-based power transformation for
#### compositional data. In Proceedings of the 4rth Compositional Data Analysis Workshop, Girona, Spain.
#### mtsagris@yahoo.gr
################################

frechet <- function(x, a) {
  ## x contains the compositional data
  ## a is the power parameter, usually between -1 and 1
  if ( a == 0 ) {
     xa <- log(x)
     y <- xa - Rfast::rowmeans(xa)
     m1 <- exp( Rfast::colmeans(y) )
     m <- m1 / sum( m1 )  ## closed geometric mean
  }  else {
     xa <- x^a
     z <- xa / Rfast::rowmeans(xa)
     m1 <- Rfast::colmeans(z) ^ ( 1 / a )
     m <- m1 / sum(m1)  ## frechet mean in general
  }
  m
}
