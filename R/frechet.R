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
  if ( length(a) == 1 ) {
    if ( a == 0 ) {
      m1 <- exp( Rfast::colmeans( Rfast::Log(x) ) )
      m <- m1 / sum( m1 )  ## closed geometric mean
    } else {
      xa <- x^a
      z <- xa / Rfast::rowmeans(xa)
      m1 <- Rfast::colmeans(z) ^ ( 1 / a )
      m <- m1 / sum(m1)  ## frechet mean in general
    }

  } else {
    m <- matrix( nrow = length(a), ncol = dim(x)[2] )
    rownames(m) <- paste("alpha=", a, sep = "")
    for ( i in 1:length(a) ) {
      if ( a[i] == 0 ) {
        m1 <- exp( Rfast::colmeans( Rfast::Log(x) ) )
        m[i, ] <- m1 / sum( m1 )  ## closed geometric mean
      } else {
        xa <- x^a[i]
        z <- xa / Rfast::rowmeans(xa)
        m1 <- Rfast::colmeans(z) ^ ( 1 / a[i] )
        m[i, ] <- m1 / sum(m1)  ## frechet mean in general
      }
    }
  }
  m
}
