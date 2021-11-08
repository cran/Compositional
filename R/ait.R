ait <- function(x, a, h = TRUE) {

  x <- as.matrix(x)
  D <- dim(x)[2]
  if (D == 1)  x <- t(x)
  D <- dim(x)[2] ## number of components

  if (a != 0) {
    z <- x^a
    z <- 1/a * ( z - Rfast::rowmeans(z) )
  } else {
    z <- Rfast::Log(x)
    z <- z - Rfast::rowmeans(z)
  }

  if (h)  z <- tcrossprod(z, helm(D))
  z
}

