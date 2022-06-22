zad.est <- function(y) {
  ## y is the compositional data
  ## x is the independent variable(s)
  dm <- dim(y)
  D <- dm[2]   ;   d <- D - 1
  ## d is the dimensionality of the simplex
  n <- dm[1] ## sample size
  x <- matrix(1, nrow = n, ncol = 1)
  runtime <- proc.time()
  ## next we separate the compositional vectors, those which contain
  ## zeros and those without. The same separation is performed for the
  ## independent variable(s)
  a1 <- which( Rfast::rowsums( y > 0 ) == D )
  a2 <- which( Rfast::rowsums( y > 0 ) != D )
  n1 <- length(a1)
  n2 <- n - n1
  ## n1 is the sample size of the compositional vectors with no zeros
  ## n2 is the sample size of the compositional vectors with zeros
  za <- y[a2, , drop = FALSE]
  za[za == 0] <- 1
  za[ za < 1 ] <- 0
  theta <- table( apply(za, 1, paste, collapse = ",") )
  theta <- as.vector(theta)
  con <- n1 * log(n1/n) + sum( theta * log(theta/n) )

  y1 <- y[a1, , drop = FALSE]
  ly1 <- log( y1 )
  x1 <- x[a1, , drop = FALSE]
  ly2 <- log( y[a2, , drop = FALSE] )
  x2 <- x[a2, , drop = FALSE]
  n1 <- nrow(y1)    ;    n2 <- n - n1
  beta.ini <- .lm.fit(x1, ly1[, -1] - ly1[, 1])$coefficients
  ini.phi <- sum( Compositional::diri.nr(y1, type = 2)$param )
  ##############
  ini.par <- c( log(ini.phi), as.vector( t( beta.ini) ) )  ## initial parameter values
  z <- list(ly1 = ly1, ly2 = ly2, x1 = x1, x2 = x2, a1 = a1, a2 = a2)
  #suppressWarnings()
  qa <- nlm( mixreg, ini.par, z = z )
  el1 <- -qa$minimum
  qa <- nlm( mixreg, qa$estimate, z = z )
  el2 <-  - qa$minimum
  qa <- optim( qa$estimate, mixreg, z = z, hessian = TRUE, control = list(maxiters = 10000) )

  phi <- exp( qa$par[1] )  ## final phi value
  mu <- c(1, exp(qa$par[-1]) ) ## final beta values
  mu <- mu / sum(mu)

  runtime <- proc.time() - runtime

  list(loglik = -qa$value + con, phi = phi, mu = mu, runtime = runtime)
}






