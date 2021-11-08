## Predictors are compositional data
alfa.knn.reg <- function(xnew, y, x, a = 1, k = 2:10, apostasi = "euclidean", method = "average") {
  nu <- dim(xnew)[1]
  
  if ( !is.null(a) ) {
    znew <- Compositional::alfa(xnew, a, h = FALSE)$aff
    z <- Compositional::alfa(x, a, h = FALSE)$aff
  } else {
    znew <- xnew
    z <- x
  }
  
  g <- Rfast::knn(xnew = znew, y = y, x = z, k = k, dist.type = apostasi, type = "R", method = method)
  colnames(g) <- paste("k=", k, sep = "")
  g
}
