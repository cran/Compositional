probout <- function(mu, su, a) {
  n <- 1e+7 
  z <- Rfast::rmvnorm(n, mu, su)
  D <- dim(z)[2]
  z1 <- z %*% helm(D + 1)
  mi <- abs(a) * Rfast::rowMins(z1, value = TRUE)
  p1 <- mean( mi <  - 1)
  z <- Rfast::rmvnorm(n, mu, su)
  z1 <- z %*% helm(D + 1)
  mi <- abs(a) * Rfast::rowMins(z1, value = TRUE)
  p2 <- mean( mi <  - 1)
  0.5 * (p1 + p2)
}
 
