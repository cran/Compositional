kl.alfapcr <- function(y, x, covar = NULL, a, k, xnew = NULL,
                       B = 1, ncores = 1, tol = 1e-07, maxiters = 50) {
  z <- Compositional::alfa(x, a)$aff
  p <- dim(z)[2]
  if (k > p)   k <- p
  eig <- prcomp(z, center = FALSE, scale = FALSE)
  vec <- eig$rotation[, 1:k, drop = FALSE]
  sc <- eig$x[, 1:k, drop = FALSE]
  if ( !is.null(covar) )  sc <- cbind(sc, covar)
  if ( !is.null(xnew) )  {
    xnew <- Compositional::alfa(xnew, a)$aff
    xnew <- cbind(xnew %*% vec, covar)
  }
  Compositional::kl.compreg(y, sc, xnew = xnew, B = B, ncores = ncores,
                            tol = tol, maxiters = maxiters)
}
