alfa.nb <- function(xnew, x, ina, a, type = "gaussian"){

  if ( !is.matrix(xnew) )  xnew <- t( as.matrix(xnew) )
  nu <- dim(xnew)[1]
  la <- length(a)
  mat <- matrix(nrow = nu, ncol = la)
  colnames(mat) <- paste("a=", a, sep = "")
  if (la > 1) {
    if ( min(x) == 0  |  min(xnew) == 0 )  a <- a[a > 0]
  }

  if ( type == "gaussian" ) {
    nb <- Rfast::gaussian.nb
  } else if ( type == "cauchy" ) {
    nb <- Rfast2::cauchy.nb
  } else if ( type == "laplace" ) {
    nb <- Rfast2::laplace.nb
  }

  for (i in 1:la) {
    y <- Compositional::alfa(x, a[i])$aff
    ynew <- Compositional::alfa(xnew, a[i])$aff
    mat[, i] <- nb(ynew, y, ina)$est
  }

  mat
}
