rgendiri <- function(n, a, b) {
  ##  n is the sample size
  ##  a is the parameters vector
  D <- length(a)
  y <- matrix( rgamma(n * D, a, b), ncol = D, byrow = TRUE )
  y / Rfast::rowsums(y)  ## Dirichlet simulated values
}
