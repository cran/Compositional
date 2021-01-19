fold.contour <- function(m, s, p, a, n = 100, x = NULL) {

  x1 <- seq( 0.001, 0.999, length = n )
  x2 <- seq( 0.001, sqrt(3)/2 - 0.001, length = n )
  mat <- matrix(NA, n, n)
  down <- ( (2 * pi) ^ (-2) ) * ( det(s) ^ (-0.5) )
  ts <- solve(s)   ;   D <- 3   ;   d <- 2
  h <- t( helm(D) )
  sqrt3 <- sqrt(3)

  for ( i in 1:c(n/2) ) {
    for ( j in 1:n ) {

      if ( x2[j] < sqrt3 * x1[i] ) {  ## This checks whether the point will lie inside the triangle
        ## The next 4 lines calculate the composition
        w3 <- 2 * x2[j] / sqrt3
        w2 <- x1[i] - x2[j] / sqrt3
        w1 <- 1 - w2 - w3
        w <- abs( c(w1, w2, w3) )
        z1 <- D / a * (w^a) / sum( w^a) - 1/a
        y1 <- z1 %*% h
        lam <- min( a * z1 ) ^ (-2)
        y2 <- as.vector( lam * y1 )
        f1 <- p * down * exp( -0.5 * c( y1 - m ) %*% ts %*% t ( t( c( y1 - m ) ) ) )
        f2 <- (1 - p) * down * lam^d * exp( -0.5 * c( y2 - m ) %*% ts %*% t( t( c( y2 - m ) ) ) )
        can <- f1 + f2

        if ( abs(can) < Inf )  mat[i, j] = can
      }

    }
  }

  for ( i in c(n/2 + 1):n ) {
    for ( j in  1:n ) {

      if ( x2[j] < sqrt3 - sqrt3 * x1[i] ) {  ## This checks whether the point will lie inside the triangle
        ## The next 4 lines calculate the composition
        w3 <- 2 * x2[j] / sqrt3
        w2 <- x1[i] - x2[j] / sqrt3
        w1 <- 1 - w2 - w3
        w <- abs( c(w1, w2, w3) )
        z1 <- D / a * (w^a) / sum(w^a) - 1/a
        y1 <- z1 %*% h
        lam <- min( a * z1 ) ^ (-2)
        y2 <- as.vector( lam * y1 )
        f1 <- p * down * exp( -0.5 * c( y1 - m ) %*% ts %*% t( t( c( y1 - m ) ) ) )
        f2 <- (1 - p) * down * lam^d * exp( -0.5 * c( y2 - m ) %*% ts %*% t( t( c( y2 - m ) ) ) )
        can <- f1 + f2

        if ( abs(can) < Inf )  mat[i, j] = can

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
