kl.diri.normal <- function(a) {
 
  d <- length(a) - 1
  da <- sum( lgamma(a) ) - lgamma( sum(a) )
  A <- sum(a)
  0.5 * d * ( 1 + log(2 * pi) ) - exp(da) + 
  sum( digamma(a) ) - A * digamma(A) + 0.5 * log( trigamma(a) ) + 
  0.5 * log( sum(1 / trigamma(a)) )

}

