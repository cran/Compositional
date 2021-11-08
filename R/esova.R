## computes ESOV distance between xnew and x
esova <- function(xnew, x) {
  dm <- dim(x)
  n <- dm[1]    ;    p <- dm[2]
  xnew <- as.matrix(xnew)
  xnew <- matrix(xnew, ncol = p ) ## makes sure xnew is a matrix
  nu <- dim(xnew)[1]
  disa <- matrix(0, n, nu)
  tx <- t(x)
  
  for (i in 1:nu) {
    xan <- xnew[i, ]
    ma <- 0.5 * ( tx + xan )
    disa[, i] <- colSums( xan * log( xan / ma ) + tx * log( tx/ma ), na.rm = TRUE )
  }
  
  disa
}
