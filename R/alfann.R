################################
#### alpha-distance
#### Tsagris Michail 5/2013
#### References: Tsagris, M. T., Preston, S., and Wood, A. T. A. (2011).
#### A data-based power transformation for
#### compositional data. In Proceedings of the 4rth Compositional Data Analysis Workshop, Girona, Spain.
#### mtsagris@yahoo.gr
################################
alfann <- function(xnew, x, a, k = 10, rann = FALSE) {
  ## x contains the compositional data
  ## a is the power parameter, usually between -1 and 1
  D <- dim(x)[2]
  if (a != 0 ) {
    y <- x^a
    ta <- Rfast::rowsums(y)
    y <- D / a * y / ta - 1/a
    ynew <- xnew^a
    ta <- Rfast::rowsums(ynew)
    ynew <- D / a * ynew / ta - 1/a
  } else {
    ya <- Rfast::Log(x)
    y <- ya - Rfast::rowmeans( ya )
    ynew <- Rfast::Log(xnew)
    ynew <- ynew - Rfast::rowmeans( ynew )
  }
  if (rann) {
    di <- RANN::nn2( data = y, query = ynew, k = k )$nn.idx
  } else  di <- Rfast::dista(ynew, y, k = k, index = TRUE, square = TRUE)
  di
}
