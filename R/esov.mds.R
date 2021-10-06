esov.mds <- function(x, k = 2, eig = TRUE) {
  d <- Compositional::esov(x)
  cmdscale(d = as.dist(d), k = k, eig = eig)
}
