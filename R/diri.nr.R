################################
#### Dirichlet distribution parameters
#### via Newton-Raphson
#### Tsagris Michail 8/2015
#### mtsagris@yahoo.gr
#### References: Estimating a Dirichlet distribution (2012)
#### Thomas P. Minka
#### http://research.microsoft.com/en-us/um/people/minka/papers/dirichlet/minka-dirichlet.pdf
################################

diri.nr <- function(x, type = 1, tol = 1e-07) {
  ## x is compositional data
  ## x can be either "lik" or "ent"

  if (type == 1) {

    runtime <- proc.time()
    x <- as.matrix(x)  ## makes sure x is a matrix
    x <- x/rowSums(x)  ## makes sure x is compositional data
    n <- nrow(x)  ## the sample size
    p <- ncol(x)  ## dimensionality
    m <- colMeans(x)
    zx <- t( log(x) )

    lm <- log(m)
    down <-  - sum(  m * rowMeans( zx - lm ) )
    sa <- 0.5 * (p - 1) / down  ## initial value for precision
    a1 <- sa * m  ## initial values
    gm <- rowSums(zx)

    z <- n * digamma( sa )
    g <- z - n * digamma(a1) + gm
    qk <-  - n * trigamma(a1)
    b <- ( sum(g / qk) ) / ( 1/z - sum(1 / qk) )
    a2 <- a1 - (g - b)/qk
    i <- 2

    while( sum( abs( a2 - a1 ) ) > tol ) {
      i <- i + 1
      a1 <- a2
      z <- n * digamma( sum(a1) )
      g <- z - n * digamma(a1) + gm
      qk <-  - n * trigamma(a1)
      b <- ( sum(g / qk) ) / ( 1/z - sum(1 / qk) )
      a2 <- a1 - (g - b) / qk
    }

    a <- a2

    loglik <- n * lgamma( sum(a) ) - n * sum( lgamma(a) ) +
      sum( zx * (a - 1) )

    runtime <- proc.time() - runtime

  } else if (type == 2) {

    runtime <- proc.time()
    x <- as.matrix(x)  ## makes sure x is a matrix
    x <- x/rowSums(x)  ## makes sure x is compositional data
    n <- nrow(x)  ## sample size
    p <- ncol(x)
    zx <- t( log(x) )

    ma <- rowMeans(zx)
    m <- colMeans(x)
    lm <- log(m)
    down <-  - sum( m * ( ma - lm ) )
    sa <- 0.5 * (p - 1) / down  ## initial value for precision
    a1 <- sa * m  ## initial values

    f <- ma - digamma(a1) + digamma( sa )
    der <-  - trigamma(a1) + trigamma( sa )
    a2 <- a1 - f / der
    i <- 2

    while ( sum( abs( a2 - a1 ) ) > tol ) {
      a1 <- a2
      i <- i + 1
      sa <- sum( a1)
      f <- ma - digamma(a1) + digamma( sa )
      der <-  - trigamma(a1) + trigamma( sa )
      a2 <- a1 - f / der

    }

    a <- a2

    loglik <- n * lgamma( sum(a) ) - n * sum( lgamma(a) ) +
      sum( zx * (a - 1) )

    runtime <- proc.time() - runtime

  }

  if ( is.null(colnames(x)) ) {
    names(a) <- paste("X", 1:p, sep = "")
  } else  names(a) <- colnames(x)

  list(iter = i, loglik = loglik, param = a, runtime = runtime)

}
