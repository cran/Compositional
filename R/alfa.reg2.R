alfa.reg2 <- function(y, x, a, xnew = NULL, yb = NULL, seb = FALSE) {
  res <- list()
  for ( i in 1:length(a) ) {
    res[[ i ]] <- Compositional::alfa.reg(y, x, a[i], xnew, yb, seb)
  }
  res
}
