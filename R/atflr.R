atflr <- function(y, x, a = seq(0.1, 1, by = 0.1), xnew ) {

  if ( min(y) == 0 )  a <- abs(a)
  res <- list()

  for ( j in 1:length(a) ) {
    ya <- y^a[j]
    ya <- ya / Rfast::rowsums(ya)
    be <- Compositional::tflr(ya, x)
    est <- ( xnew %*% be )^( 1/a[j] )
    est <- est / Rfast::rowsums(est)
    res[[ j ]] <- est
  }

  res
}

