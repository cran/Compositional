################################
#### Contour plot of the Dirichlet distribution in S^2
#### Tsagris Michail 1/2013
#### mtsagris@yahoo.gr
################################
fd.contour <- function(alpha, prob, tau, n = 100, x = NULL) {
  ## a are the estimated Dirichlet parameters
  ## n shows the number of points at which the density is calculated
  ## so, n^2 points are used.
  ## x should be a 3-part compositional data or NULL for no data
  x1 <- seq(0.001, 0.999, length = n)  ## coordinates of x
  sqrt3 <- sqrt(3)
  x2 <- seq(0.001, sqrt3/2 - 1e-03, length = n)  ## coordinates of y
  mat <- matrix(nrow = n, ncol = n)

  for ( i in 1:c(n/2) ) {
    for (j in 1:n) {
      if ( x2[j] < sqrt3 * x1[i] ) {
        ## This checks if the point will lie inside the triangle
        ## the next three lines invert the points which lie inside
        ## the triangle back into the composition in S^2
        w3 <- 2 * x2[j] / sqrt3
        w2 <- x1[i] - x2[j]/sqrt3
        w1 <- 1 - w2 - w3
        w <- c(w1, w2, w3)
        can <- Compositional::fd.density(w, alpha, prob, tau)
        if (abs(can) < Inf)  mat[i, j] <- can
      }
    }
  }
  for (i in c(n/2 + 1):n) {
    for (j in 1:n) {
      ## This checks if the point will lie inside the triangle
      if ( x2[j] < sqrt3 - sqrt3 * x1[i] ) {
        ## the next three lines invert the points which lie inside
        ## the triangle back into the composition in S^2
        w3 <- 2 * x2[j] / sqrt3
        w2 <- x1[i] - x2[j]/sqrt3
        w1 <- 1 - w2 - w3
        w <- c(w1, w2, w3)
        can <- Compositional::fd.density(w, alpha, prob, tau)
        if (abs(can) < Inf)  mat[i, j] <- can
      }
    }
  }
  contour(x1, x2, mat, col = 3, xlab = " ", ylab = " ",
          pty = "s", xaxt = "n", yaxt = "n", bty = "n")
  b1 <- c(0.5, 0, 1, 0.5)
  b2 <- c(sqrt3/2, 0, 0, sqrt3/2)
  b <- cbind(b1, b2)
  ## the next line draws the triangle in the two dimensions
  points(b[, 1], b[, 2], type = "l", xlab = " ", ylab = " ")

  if ( !is.null(x) ) {
    proj <- matrix(c(0, 1, 0.5, 0, 0, sqrt3/2), ncol = 2)
    xa <- x %*% proj
    points(xa[, 1], xa[, 2])
    nam <- colnames(x)
    if ( is.null(nam) )  nam <- paste("X", 1:3, sep = "")
    points(b[, 1], b[, 2], type = "l", xlab = " ", ylab = " ")
    text( b[1, 1], b[1, 2] + 0.02, nam[3], cex = 1 )
    text( b[2:3, 1], b[2:3, 2] - 0.02, nam[1:2], cex = 1 )
  }

}

