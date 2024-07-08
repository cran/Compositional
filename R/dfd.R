dfd <- function(x, alpha, prob, tau) {
  if ( is.null(dim(x)[1]) )  x <- matrix(x, nrow = 1)
  aplus <- sum(alpha)
  f1 <- lgamma(aplus + tau) - sum( lgamma(alpha) ) + as.vector( log(x) %*% (alpha - 1) ) + tau * log(x)
  f2 <- log(prob) + lgamma(alpha) - lgamma(alpha + tau)
  Rfast::rowsums( exp( Rfast::eachrow(f1, f2, oper = "+") ) )
}


#dfd <- function(x, alpha, prob, tau){
#  FlexDir::FD.density(x = x, a = alpha, p = prob, t = tau)
#}
