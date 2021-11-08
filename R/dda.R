dda <- function(xnew, x, ina) {
  if ( !is.matrix(xnew) )  xnew <- matrix(xnew, nrow = 1)
  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  mat <- matrix(nrow = dim(xnew)[1], ncol = g)
  
  for (j in 1:g) {
    a <- Compositional::diri.nr(x[ina == j, ])$param
    mat[, j] <- Compositional::ddiri(xnew, a, logged = TRUE)
  }
  
  Rfast::rowMaxs(mat)
}
