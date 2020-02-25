divergence <- function(x, type = "kullback_leibler", vector = FALSE) {
  Rfast::Dist(x, method = type, vector = vector)
}
