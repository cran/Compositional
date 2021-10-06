alfa.mds <- function(x, a, k = 2, eig = TRUE) {
  d <- Compositional::alfadist(x, a = a)
  cmdscale(d = as.dist(d), k = k, eig = eig)
}
