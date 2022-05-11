a.mle <- function(a, x) {
  ## x is the compositional data
  ## a is the value of the alpha parameter
  dm <- dim(x)
  n <- dm[1]    ;    D <- dm[2]   ## dimensions of x
  d <- D - 1  ## dimensionality of the simplex
  ja <- sum( Rfast::Log(x) )  ## part of the Jacobian determinant
  #########
  if ( abs(a) < 1e-9 ) {
    aff <- Compositional::alef(x, 0)$aff
    su <- Rfast::cova(aff)
    loglik <-  - 0.5 * n * d - 0.5 * n * log( abs( det(2 * pi * (n - 1)/n * su) ) ) - ja - n/2 * log(D)

  } else {
    mod <- Compositional::alef(x, a)
    y <- mod$aff
    sk <- mod$sk
    lam <- 1 / ( a^2 * Rfast::rowMins(y, value = TRUE)^2 )    ##  1/apply(a * y, 1, min)^2
    y1 <- y %*% t( helm(D) )
    y2 <- y1 * lam
    com <-  - 0.5 * n * d * log(2 * pi ) + n * (d + 0.5) * log(D) + (a - 1) * ja - D * sum( log(sk) ) + 0.5 * d
    ## step 1
    ma <- Rfast::colmeans(y1)
    sa <- Rfast::cova(y1)
    f1 <- exp( -0.5 * Rfast::mahala(y1, ma, sa) )
    f2 <- lam^d * exp(-0.5 * Rfast::mahala(y2, ma, sa) )
    p <- f1/(f1 + f2)
    per <- sum(p) / n
    ela1 <- sum( log(per * f1 + (1 - per) * f2), na.rm = TRUE ) - 0.5 * n * log( det(sa) )
    ma <- colMeans(p * y1 + (1 - p) * y2, na.rm = TRUE)
    z1 <- sqrt(p) * ( y1 - rep(ma, rep(n, d)) )
    z1 <- na.omit(z1)
    z2 <- sqrt(1 - p) * ( y2 - rep(ma, rep(n, d)) )
    z2 <- na.omit(z2)
    sa <- ( crossprod(z1) + crossprod(z2) )/n
    f1 <- exp( -0.5 * Rfast::mahala(y1, ma, sa) )
    f2 <- lam^d * exp( -0.5 * Rfast::mahala(y2, ma, sa) )
    p <- f1/(f1 + f2)
    per <- sum(p) / n
    ela2 <- sum( log(per * f1 + (1 - per) * f2), na.rm = TRUE ) - 0.5 * n * log( det(sa) )

    while ( abs(ela2 - ela1 ) > 1e-06 ) {
      ela1 <- ela2
      ma <- colMeans(p * y1 + (1 - p) * y2, na.rm = TRUE)
      z1 <- sqrt(p) * ( y1 - rep(ma, rep(n, d)) )
      z1 <- na.omit(z1)
      z2 <- sqrt(1 - p) * ( y2 - rep(ma, rep(n, d)) )
      z2 <- na.omit(z2)
      sa <- ( crossprod(z1) + crossprod(z2) )/n
      f1 <- exp( -0.5 * Rfast::mahala(y1, ma, sa) )
      f2 <- lam^d * exp( -0.5 * Rfast::mahala(y2, ma, sa) )
      p <- f1 / (f1 + f2)
      per <- sum(p) / n
      ela2 <- sum( log(per * f1 + (1 - per) * f2), na.rm = TRUE ) - 0.5 * n * log( det(sa) )
    }  ## end while ( abs(ela2 - ela1 ) > 1e-06 )
    loglik <- ela2 + com
  }  ## end if ( abs(a) < 1e-9)

  loglik
}




