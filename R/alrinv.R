alrinv <- function(y) {
  x <- cbind(1, exp(y) )
  x / rowsums(x)
}
