mixreg <- function(param, z) {
  ## separation of phi and the betas
  ## and the exponential to avoid negative values of phi
  phi <- exp( param[1] )   ;   para <- param[-1]
  ly1 <- z$ly1     ;     ly2 <- z$ly2
  x1 <- z$x1     ;     x2 <- z$x2
  a1 <- z$a1     ;     a2 <- z$a2
  dm <- dim(ly1)
  D <- dm[2]    ;    d <- D - 1
  n1 <- length(a1)    ;   n2 <- length(a2)
  ## n1 is the sample size of the compositional vectors with no zeros
  ## n2 is the sample size of the compositional vectors with zeros
  n <- n1 + n2  ## total sample size
  ## next we separate the compositional vectors, those which contain
  ## zeros and those without. The same separation is performed for the
  ## independent variable(s)
  be <- matrix(para, ncol = d)   ## be is the matrix of the betas
  be <- cbind(0, be)
  mu1 <- exp(x1 %*% be)
  mu <- mu1 / Rfast::rowsums(mu1) ## fitted values
  ## next we find the fitted values for the compositional vectors with zeros
  ly3 <- ly2
  ind <- which(is.infinite(ly2))
  ly3[ind] = 0
  mu2 <- exp(x2 %*% be )
  mu2[ind] <- 0
  mu2 <- mu2/rowsums(mu2)
  zeros <- - sum( lgamma(phi * mu2[-ind]) ) + sum( (mu2 * phi - 1) * ly3 )
  ba <- phi * mu
  f <-  - n * lgamma(phi) + sum( lgamma( ba[ba>0] ) ) - sum( (ba - 1) * ly1 ) - zeros
  f
}




