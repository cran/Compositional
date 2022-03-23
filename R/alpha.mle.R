alpha.mle <- function(x, a) {
  ## x is the compositional data
  ## a is the value of the alpha parameter
  dm <- dim(x)
  n <- dm[1]  ;  D <- dm[2]  ## dimensions of x
  d <- D - 1  ## dimensionality of the simplex
  ja <- sum( Rfast::Log(x) )  ## part of the Jacobian determinant
  #########
  if ( abs(a) < 1e-9 ) {  ## i.e. if alpha = 0
    mod <- alfa(x, 0)
    aff <- mod$aff
    su <- Rfast::cova(aff)
    con <-  - n/2 * d * log(2 * pi * (n - 1)/n ) - (n - 1) * d/2 + n * (d + 0.5) * log(D)
    lik <-  - n/2 * log( abs( det( cov(aff) ) ) ) - ja - D * mod$sa + con
    result <- list(loglik = lik, mu = Rfast::colmeans(aff), su = su)

  } else {
    mod <- alef(x, a)
    y <- mod$aff
    sk <- mod$sk
    lam <- 1 /(a^2 * Rfast::rowMins(y, value = TRUE)^2)   ##  1 / apply(a * y, 1, min)^2
    y1 <- y %*% t( helm(D) )
    y2 <- y1 * lam
    lamd <- lam^d
    ma <- Rfast::colmeans(y1)
    sa <- Rfast::cova(y1)
    com <-  - 0.5 * n * d * log(2 * pi ) + n * (d + 0.5) * log(D) + (a - 1) * ja - D * sum( log(sk) )
    ## step 1
    con <-  - 0.5 * n * log( det(sa) )
    f1 <- exp( -0.5 * Rfast::mahala(y1, ma, sa) )
    f2 <- lamd * exp(-0.5 * Rfast::mahala(y2, ma, sa) )
    p <- f1 / (f1 + f2)
    per <- sum(p) / n
    ela1 <- sum( log(per * f1 + (1 - per) * f2) ) + con
    ## step 2
    ma <- colMeans(p * y1 + (1 - p) * y2, na.rm = TRUE)
    z1 <- sqrt(p) * Rfast::eachrow(y1, ma, oper = "-")       ## ( y1 - rep(ma, rep(n, d)) )
	  z2 <- sqrt(1 - p) * Rfast::eachrow(y2, ma, oper = "-")   ## ( y2 - rep( ma, rep(n, d) ) )
    sa <- ( crossprod(z1) + crossprod(z2) )/n
    con <-  - 0.5 * n * log( det(sa) )
    f1 <- exp( -0.5 * Rfast::mahala(y1, ma, sa) )
    f2 <- lamd * exp( -0.5 * Rfast::mahala(y2, ma, sa) )
    p <- f1 / (f1 + f2)
    per <- sum(p) / n
    ela2 <- sum( log(per * f1 + (1 - per) * f2) ) + con
    ## step 3 and beyond
    k <- 2
    while ( abs(ela2 - ela1) > 1e-06 ) {
     k <- k + 1
     ela1 <- ela2
     ma <- colMeans(p * y1 + (1 - p) * y2, na.rm = TRUE)
     z1 <- sqrt(p) * Rfast::eachrow(y1, ma, oper = "-")       ## ( y1 - rep(ma, rep(n, d)) )
	   z2 <- sqrt(1 - p) * Rfast::eachrow(y2, ma, oper = "-")   ## ( y2 - rep( ma, rep(n, d) ) )
     sa <- ( crossprod(z1) + crossprod(z2) )/n
     con <-  - 0.5 * n * log( det(sa) )
     f1 <- exp( -0.5 * Rfast::mahala(y1, ma, sa) )
     f2 <- lamd * exp(-0.5 * Rfast::mahala(y2, ma, sa) )
     p <- f1 / (f1 + f2)
     per <- sum(p) / n
     ela2 <- sum( log(per * f1 + (1 - per) * f2) ) + con
    }
   result <- list(iters = k, p = per, loglik = ela2 + com + 0.5 * d, mu = ma, su = sa)
   }
   result
}
