################################
#### alpha-distance
#### Tsagris Michail 5/2013
#### References: Tsagris, M. T., Preston, S., and Wood, A. T. A. (2011).
#### A data-based power transformation for
#### compositional data. In Proceedings of the 4rth Compositional Data Analysis Workshop, Girona, Spain.
#### mtsagris@yahoo.gr
################################
alfadista <- function(xnew, x, a, type = "euclidean", square = FALSE) {
  ## x contains the compositional data
  ## a is the power parameter, usually between -1 and 1
  y <- Compositional::alfa(x, a, h = FALSE)$aff
  ynew <- Compositional::alfa(xnew, a, h = FALSE)$aff
  Rfast::dista(ynew, y, type = type, square = square)
}
