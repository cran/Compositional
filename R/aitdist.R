aitdist <- function (x, a, type = "euclidean", square = FALSE) {
  y <- Compositional::ait(x, a, h = FALSE)
  Rfast::Dist(y, method = type, square = square)
}
