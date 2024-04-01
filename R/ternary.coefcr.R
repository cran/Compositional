ternary.coefcr <- function(y, x, type = "scls", conf = 0.95, R = 1000, dg = FALSE, hg = FALSE) {

  dm <- dim(x)
  n <- dm[1]  ;  px <- dm[2]
  bb <- matrix( nrow = R, ncol = 3 * px )

  if ( type == "scls" ) {
    B <- Compositional::scls(y, x)$be
    for (i in 1:R) {
      id <- Rfast2::Sample.int(n, n, replace = TRUE)
      bb[i, ] <- as.vector( t( Compositional::scls(y[id, ], x[id, ])$be ) )
    }
  } else {
    B <- Compositional::tflr(y, x)$be
    for (i in 1:R) {
      id <- Rfast2::Sample.int(n, n, replace = TRUE)
      bb[i, ] <- as.vector( t( Compositional::tflr(y[id, ], x[id, ])$be ) )
    }
  }

  bb <- abs(bb)
  confr <- NULL
  ind <- matrix(1:(3 * px), nrow = 3)
  for (i in 1:px) {
    fit <- MASS::cov.rob(bb[, ind[1:2, i]], quantile.used = ceiling(conf * R), method = "mve")
    cr <- predict( cluster::ellipsoidhull( bb[fit$best, ind[1:2, i] ] ) )
    cr[cr < 0] <- 0
    cr[cr > 1] <- 1
    cr3 <- pmax(0, 1 - Rfast::rowsums(cr) )
    cr <- cbind(cr, cr3)
    cr <- cr / Rfast::rowsums(cr)
    confr <- cbind(confr, cr)
  }

  n <- dim(confr)[1]
  ina <- numeric(n) + 20
  ## the next code checks for zeros
  b1 <- c(0.5, 0, 1, 0.5)
  b2 <- c(sqrt(3)/2, 0, 0, sqrt(3)/2)
  b <- cbind(b1, b2)
  plot(b[, 1], b[, 2], type = "l", xlab = " ", ylab = " ", pty = "s",
       xaxt = "n", yaxt = "n", bty = "n", lwd = 2)
  proj <- matrix(c(0, 1, 0.5, 0, 0, sqrt(3)/2), ncol = 2)

  text( b[1, 1], b[1, 2] + 0.02, "Y3", col = "black", font = 2 )
  text( b[2, 1] + 0.02, b[2, 2] - 0.02, "Y1", col = "black", font = 2 )
  text( b[3, 1] - 0.02, b[2, 2] - 0.02, "Y2", col = "black", font = 2 )

  for (i in 1:px) {
    d <- confr[, ind[, i]] %*% proj
    lines( d[1:n, 1], d[1:n, 2], lwd = 2, col = i)
  }
  d <- B %*% proj
  points(d, pch = 3, lwd = 2, col = 1:px)

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

}
