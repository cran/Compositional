alfa.pca <- function(x, a, center = TRUE, scale = TRUE, k = NULL, vectors = FALSE) {
  x <- Compositional::alfa(x, a, h = FALSE)$aff
  Rfast2::pca(x, center, scale, k, vectors )
}
