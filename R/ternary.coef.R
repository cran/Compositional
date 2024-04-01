ternary.coef <- function(B, dg = FALSE, hg = FALSE, colour = NULL) {

  if ( !is.null( colnames(B) ) ) {
    nam <- colnames(B)
  } else nam <- paste("Y", 1:3, sep = "")

  n <- dim(B)[1]
  b1 <- c(0.5, 0, 1, 0.5)
  b2 <- c(sqrt(3)/2, 0, 0, sqrt(3)/2)
  b <- cbind(b1, b2)
  plot(b[, 1], b[, 2], type = "l", xlab = " ", ylab = " ", pty = "s",
  xaxt = "n", yaxt = "n", bty = "n", lwd = 2)
  proj <- matrix(c(0, 1, 0.5, 0, 0, sqrt(3)/2), ncol = 2)
  d <- B %*% proj
  if ( is.null(colour) )  colour <- numeric(n) + 1

  text( b[1, 1], b[1, 2] + 0.02, nam[3], col = "black", font = 2 )
  text( b[2, 1] + 0.02, b[2, 2] - 0.02, nam[1], col = "black", font = 2 )
  text( b[3, 1] - 0.02, b[2, 2] - 0.02, nam[2], col = "black", font = 2 )

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

  for (i in 1:n)  text(d[i, 1], d[i, 2],  substitute( B[i], list(i = i) ), col = colour)

}


