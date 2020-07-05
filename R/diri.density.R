diri.density <- function(x, a, logged = FALSE) {
  if ( !is.matrix(x) )  {
    x <- as.matrix(x)
    x <- t(x)
  }
  f <- lgamma( sum(a) ) - sum( lgamma(a) ) + as.numeric( log(x) %*% (a - 1) )
  if ( logged )  f <- exp(f)
  f
}
