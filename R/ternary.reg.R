ternary.reg <- function(y, est, id, labs) {

  if ( !is.null( colnames(y) ) ) {
    nam <- colnames(y)
  } else nam <- paste("y", 1:3, sep = " ")

  n <- dim(y)[1]
  idd <- as.numeric(id)
  k <- max(id)
  ## the next code checks for zeros
  ina <- rep(1, n)
  ina[ rowSums(y == 0) > 0 ] <- 3
  b1 <- c(0.5, 0, 1, 0.5)
  b2 <- c(sqrt(3)/2, 0, 0, sqrt(3)/2)
  b <- cbind(b1, b2)
  plot(b[, 1], b[, 2], type = "l", xlab = " ", ylab = " ", pty = "s",
       xaxt = "n", yaxt = "n", bty = "n", lwd = 2)
  proj <- matrix(c(0, 1, 0.5, 0, 0, sqrt(3)/2), ncol = 2)
  d <- y %*% proj
  points( d[, 1], d[, 2], col = ina )
  text( b[1, 1], b[1, 2] + 0.02, nam[3], col = "black", font = 2 )
  text( b[2:3, 1], b[2:3, 2] - 0.02, nam[1:2], col = "black", font = 2)

  d <- est %*% proj
  for (i in 1:k)  {
    d1 <- d[idd == i, ]
    d1 <- d1[ order(d1[, 1]), ]
    points(d1[, 1], d1[, 2], pch = 16, col = i)
  }

  legend("topright", labs, lwd = rep(2, k), col = unique(id), text.col = unique(id), bg = 'gray90')

}
