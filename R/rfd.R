rfd <- function(n, alpha, prob, tau) {
  D <- length(alpha)
  x <- t( rmultinom(n, 1, prob) )
  ind <- Rfast::rowMaxs(x)
  x <- x * rgamma(n, alpha[ind] + tau)

  for (j in 1:D) {
   id <- which( x[, j] == 0 )
   x[id, j] <- rgamma(length(id), alpha[j])
  }
  x / Rfast::rowsums(x)
}

#rfd <- function(n, alpha, prob, tau) {
#  FlexDir::FD.generate(n = n, a = alpha, p = prob, t = tau)
#}
