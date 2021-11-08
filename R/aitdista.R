aitdista <- function(xnew, x, a, type = "euclidean", square = FALSE) {
  ## x contains the compositional data
  ## a is the power parameter, usually between -1 and 1
  y <- Compositional::ait(x, a, h = FALSE)
  ynew <- Compositional::ait(xnew, a, h = FALSE)
  Rfast::dista(ynew, y, type = type, square = square)
}
