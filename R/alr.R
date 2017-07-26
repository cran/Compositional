alr <- function(x) {
  Rfast::Log( x[, -1]/x[, 1] ) 
}
