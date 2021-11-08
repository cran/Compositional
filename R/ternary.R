################################
#### Ternary plot
#### Tsagris Michail 5/2012
#### mtsagris@yahoo.gr
################################

ternary <- function(x, dg = FALSE, hg = FALSE, means = TRUE, pca = FALSE, colour = NULL) {
  ## x contains the compositional data
  ## if means==TRUE it will plot the arithmetic and the
  ## closed geometric mean
  ## if pca==TRUE it will plot the first principal component
  if ( !is.null( colnames(x) ) ) {
    nam <- colnames(x)
  } else nam <- paste("X", 1:3, sep = "")

  n <- dim(x)[1]
  ina <- numeric(n) + 20
  ## m1 is the closed geometric mean
  g1 <- Rfast::colmeans( log(x[, -1] / x[, 1]) )
  g2 <- c( 1, exp(g1) )
  m1 <- g2 / sum(g2)
  ## m2 is the simple arithmetic mean
  m2 <- Rfast::colmeans(x)
  x <- rbind(x, m1, m2)
  ## the next code checks for zeros
  b1 <- c(0.5, 0, 1, 0.5)
  b2 <- c(sqrt(3)/2, 0, 0, sqrt(3)/2)
  b <- cbind(b1, b2)
  plot(b[, 1], b[, 2], type = "l", xlab = " ", ylab = " ", pty = "s",
  xaxt = "n", yaxt = "n", bty = "n")
  proj <- matrix(c(0, 1, 0.5, 0, 0, sqrt(3)/2), ncol = 2)
  d <- x %*% proj
  if ( is.null(colour) )  colour <- numeric(n) + 1
  points( d[1:n, 1], d[1:n, 2], col = colour )
  text( b[1, 1], b[1, 2] + 0.02, nam[3], cex = 1 )
  text( b[2:3, 1], b[2:3, 2] - 0.02, nam[1:2], cex = 1 )

  if ( means ) {
    ## should the means appear in the plot?
    points( d[c(n + 1), 1], d[c(n + 1), 2], pch = 2, col = 2, lwd = 2 )
    points( d[c(n + 2), 1], d[c(n + 2), 2], pch = 3, col = 3, lwd = 2 )
    legend("topright", c("closed geometric mean"," arithmetic mean"),
    pch = c(2, 3), col = c(2, 3), bg = 'gray90')
  }

  if (pca  &  min(x) > 0 ) {
    ## should the first principal component appear?
    zx <- log(x[1:n, ])
    z <- zx - Rfast::rowmeans( zx )  ## clr transformation
    m <- Rfast::colmeans(z)  ## mean vector in the clr space
    a <- eigen( Rfast::cova(z) )$vectors[, 1] + m  ## move the unit vector a bit
    sc <- z %*% a
    lam <- seq( min(sc) - 1.5, max(sc) + 1.5, length = n )
    x1 <- cbind( a[1] * lam, a[2] * lam, a[3] * lam) + cbind( m[1] * (1 - lam),
    m[2] * (1 - lam), m[3] * (1 - lam) )
    expx1 <- exp(x1)
    wa1 <- expx1 / Rfast::rowsums( expx1 )  ## first principal component in S^2
    wa <- wa1 %*% proj
    lines(wa, lwd = 2, lty = 2)
  }

  if ( dg ) {
    a1 <- matrix(0, nrow = 11, ncol = 3)
    a1[, 2] <- seq(0, 1, by = 0.1)
    a1[, 3] <- seq(1, 0, by = -0.1)
    ## right
    b1 <- a1 %*% proj
    ## left
    a2 <- cbind(a1[, 2], a1[, 1], a1[, 3])
    b2 <- a2 %*% proj
    ## horizontal
    a3 <- cbind(a1[, 2], a1[, 3], a1[, 1])
    b3 <- a3 %*% proj
    
    for ( i in 2:dim(b1)[1] ) {
      segments(x0 = b1[i, 1], y0 = b1[i, 2], x1 = b3[12 - i, 1], y1 = b3[i, 2], col = "lightgrey", lty = 2)
    }
    for (i in 1:(dim(b1)[1] - 1 ) ) {
      segments(x0 = b2[i, 1], y0 = b2[i, 2], x1 = b3[i, 1], y1 = b3[12 - i, 2], col = "lightgrey", lty = 2)
    }
    lines(b[, 1], b[, 2])
    for ( i in 2:( dim(b1)[1] - 1 ) )  {
      text(b1[i, 1] + 0.025, b1[i, 2] + 0.025, a1[i, 3], cex = 1)
      text(b2[i, 1] - 0.025, b2[i, 2] + 0.02, a1[i, 2], cex = 1)
      text(b3[i, 1], b3[i, 2] - 0.02, a1[i, 3], cex = 1)
    }
  } ## end if dg

  if ( hg ) {
    a1 <- matrix(0, nrow = 11, ncol = 3)
    a1[, 2] <- seq(0, 1, by = 0.1)
    a1[, 3] <- seq(1, 0, by = -0.1)
    ## right
    b1 <- a1 %*% proj
    ## left
    a2 <- cbind(a1[, 2], a1[, 1], a1[, 3])
    b2 <- a2 %*% proj
    for ( i in 2:c( dim(b1)[1] - 1) ) {
      segments(x0 = b1[i, 1], y0 = b1[i, 2], x1 = b2[i, 1], y1 = b2[i, 2], col = "lightgrey", lty = 2)
    }
    lines(b[, 1], b[, 2])
    for ( i in 2:( dim(b1)[1] - 1 ) )  {
      text(b1[i, 1] + 0.025, b1[i, 2] + 0.025, a1[i, 3], cex = 1)
      text(b2[i, 1] - 0.025, b1[i, 2] + 0.02, a1[i, 2], cex = 1)
    }
  } ## end if hg

  mu <- rbind(m1, m2)
  rownames(mu) <- c("closed geometric", "arithmetic mean")
  colnames(mu) <- nam
  mu
}
