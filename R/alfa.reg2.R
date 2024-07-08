alfa.reg2 <- function(y, x, a, xnew = NULL) {
  res <- list()
  for ( i in 1:length(a) ) {
    res[[ i ]] <- Compositional::alfa.reg(y, x, a[i], xnew)
  }
  res
}
