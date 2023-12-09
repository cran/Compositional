ait.knn <- function (xnew, x, ina, a = 1, k = 5, mesos = TRUE,
    apostasi = "euclidean", rann = FALSE) {

  p <- dim(x)[2]
  xnew <- as.matrix(xnew)
  xnew <- matrix(xnew, ncol = p)
  ina <- as.numeric(ina)
  nc <- max(ina)
  nu <- dim(xnew)[1]

  if (!is.null(a)) {
    znew <- Compositional::ait(xnew, a, h = FALSE)
    z <- Compositional::ait(x, a, h = FALSE)
  } else {
    znew <- xnew
    z <- x
  }

  if (rann) {
    klen <- length(k)
    di <- Rnanoflann::nn(data = x, points = xnew, k = max(k), square = TRUE)$indices
    g <- matrix(nrow = nu, ncol = klen)
    m1 <- matrix(nrow = max(k), ncol = nu)
    for (i in 1:nu)  m1[, i] <- ina[di[i, ]]
      for (j in 1:klen) g[, j] <- Rfast::colMaxs( Rfast::colTabulate(m1[1:k[j], ]) )
    } else  g <- Rfast::knn(xnew = znew, y = ina, x = z, k = k,
  dist.type = apostasi, type = "C", freq.option = 1)

  colnames(g) <- paste("k=", k, sep = "")
  g
}
