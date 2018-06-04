################################
#### Contour plot of the kernel density estimate in S^2
#### Tsagris Michail 2/2015
#### mtsagris@yahoo.gr
################################
comp.kerncontour <- function(x, type = "alr", n = 100) {
  ## x contains the compositional data
  ## type determines which log-ratio transformation will be used.
  ## If type='alr' (the default) the additive
  ## log-ratio transformation is used.
  ## If type='ilr', the isometric log-ratio is used
  ## n is the number of points of each axis used
  nu <- dim(x)[1]  ## sample size
  sqrt3 <- sqrt(3)
  ha <- t( helm(3) )

  if (type == "alr")  z <- log(x[, -3]/x[, 3])  ## alr transformation
  if (type == "ilr") {  ## isometric log-ratio transformation
      zx <- log(x)
      z <- zx - Rfast::rowmeans( zx )
      z <- z %*% ha
  }

  hopt <- mkde.tune(z)$hopt
  con <- hopt^2
  ts <- diag( 1/hopt^2, 2 )
  x1 <- seq(0.001, 0.999, length = n)
  x2 <- seq(0.001, sqrt3/2 - 0.001, length = n)
  mat <- matrix(nrow = n, ncol = n)

  for ( i in 1:c(n/2) ) {
    for ( j in 1:n ) {
      if (x2[j] < sqrt3 * x1[i]) {
        ## This checks if the point will lie inside the triangle
        ## The next 4 lines calculate the composition
        w3 <- 2 * x2[j] / sqrt3
        w2 <- x1[i] - x2[j]/sqrt3
        w1 <- 1 - w2 - w3
        w <- c(w1, w2, w3)
        if (type == "alr")  y <- log(w[-3]/w[3])  ## alr transformation
        if (type == "ilr") {  ## isometric log-ratio transformation
            y <- log(w) - mean(log(w))
            y <- as.vector(y %*% ha )
        }
        a <- numeric(nu)
        for (l in 1:nu)   a[l] <- as.vector( t(z[l, ] - y) %*% ts %*% ( z[l, ] - y ) )
        can <- 0.5 / pi / con * sum( exp(-0.5 * a) )/nu
        if ( abs(can) < Inf )  mat[i, j] <- can        
      }
    }
  }

  for ( i in c(n/2 + 1):n ) {
    for ( j in 1:n ) {
      ## This checks if the point will lie inside the triangle
      if (x2[j] < sqrt3 - sqrt3 * x1[i]) {
        ## The next 4 lines calculate the composition
        w3 <- 2 * x2[j] / sqrt3
        w2 <- x1[i] - x2[j]/sqrt3
        w1 <- 1 - w2 - w3
        w <- c(w1, w2, w3)
        if (type == "alr") y <- log(w[-3]/w[3])  ## alr transformation
        if (type == "ilr") {  ## isometric log-ratio transformation
          y <- log(w) - mean(log(w))
          y <- as.vector( y %*% ha )
        }
        a <- numeric(nu)
        for (l in 1:nu) a[l] <- as.vector( t(z[l, ] - y ) %*% ts %*% (z[l, ] - y) )
        can <- 0.5 / pi / con * sum( exp(-0.5 * a) )/nu
        if (abs(can) < Inf)  mat[i, j] <- can 
      }
    }
  }

  contour( x1, x2, mat, col = 3, pty = "s", xaxt = "n", yaxt = "n", bty = "n" )
  proj <- matrix(c(0, 1, 0.5, 0, 0, sqrt3/2), ncol = 2)
  da <- x %*% proj
  points(da[, 1], da[, 2])
  b1 <- c(0.5, 0, 1, 0.5)
  b2 <- c(sqrt3/2, 0, 0, sqrt3/2)
  b <- cbind(b1, b2)
  points(b[, 1], b[, 2], type = "l", xlab = " ", ylab = " ")
  nam <- colnames(x)
  if ( is.null(nam) )  nam <- paste("X", 1:3, sep = "")
  text( b[1, 1], b[1, 2] + 0.02, nam[3], cex = 1 )
  text( b[2:3, 1], b[2:3, 2] - 0.02, nam[1:2], cex = 1 )
}
