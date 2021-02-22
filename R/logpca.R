logpca <- function(x, center = TRUE, scale = TRUE, k = NULL, vectors = FALSE) {
  Rfast2::pca( Rfast::Log(x), center, scale, k, vectors )
}
