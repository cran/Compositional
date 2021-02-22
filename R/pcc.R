pcc <- function(x) {
  s <- cov( log(x) )
  d <- dim(s)[1]
  down <- diag(s)
  denom <- Rfast::colsums( Rfast::comb_n(down, 2) )
  denom <- Rfast::squareform(denom)
  pcc <- 2 * s / denom
  diag(pcc) <- 1
  pcc
}
