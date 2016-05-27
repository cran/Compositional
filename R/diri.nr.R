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
    ell <- NULL
    m <- colMeans(x)
    zx <- log(x)

    lm <- log(m)
    down <-  - sum(  m * colMeans( zx - rep(lm, rep(n, p ) ) ) )
    s <- 0.5 * (p - 1) / down  ## initial value for precision
    a <- s * m  ## initial values
    ell[1] <- n * lgamma( sum(a) ) - n * sum( lgamma(a) ) +
    sum( zx %*% (a- 1) )
    gm <- colSums(zx)
    z <- n * digamma( sum(a) )
    g <- z - n * digamma(a) + gm
    qk <-  - n * trigamma(a)
    b <- ( sum(g / qk) ) / ( 1/z - sum(1 / qk) )
    a <- a - (g - b)/qk
    ell[2] <- n * lgamma( sum(a) ) - n * sum( lgamma(a) ) +
    sum( zx %*% (a - 1) )
    i <- 2

    while( abs( ell[i] - ell[i - 1] ) > tol ) {
      i <- i + 1
      z <- n * digamma( sum(a) )
      g <- z - n * digamma(a) + gm
      qk <-  - n * trigamma(a)
      b <- ( sum(g / qk) ) / ( 1/z - sum(1 / qk) )
      a <- a - (g - b) / qk
      ell[i] <- n * lgamma( sum(a) ) - n * sum( lgamma(a) ) +
      sum( zx %*% (a - 1) )
    }

    loglik <- ell[i]

    runtime <- proc.time() - runtime

  } else if (type == 2) {

    runtime <- proc.time()
    x <- as.matrix(x)  ## makes sure x is a matrix
    x <- x/rowSums(x)  ## makes sure x is compositional data
    n <- nrow(x)  ## sample size
    p <- ncol(x)
    zx <- log(x)

    ma <- colMeans(zx)
    m <- colMeans(x)
    lm <- log(m)
    down <-  - sum(  m * colMeans( zx - rep(lm, rep(n, p ) ) ) )
    s <- 0.5 * (p - 1) / down  ## initial value for precision
    a <- s * m  ## initial values

    f <- ma - digamma(a) + digamma( sum(a) )
    der <-  - trigamma(a) + trigamma( sum(a) )
    der <- diag( 1 / der )
    a <- rbind(a, a - f %*% der )
    i <- 2

    while ( sum( abs( a[i, ] - a[i - 1, ] ) ) > tol ) {
      i <- i + 1
      f <- ma - digamma(a[i - 1, ]) + digamma( sum(a[i - 1, ]) )
      der <-  - trigamma(a[i - 1, ]) + trigamma( sum(a[i - 1, ]) )
      der <- diag( 1 / der )
      a <- rbind(a, a[i - 1, ] - f %*% der )

    }

    a <- a[i, ]
    loglik <- n * lgamma( sum(a) ) - n * sum( lgamma(a) ) +
    sum( zx %*% (a - 1) )

    runtime <- proc.time() - runtime

  }

  if ( is.null(colnames(x)) ) {
    names(a) <- paste("X", 1:p, sep = "")
  } else  names(a) <- colnames(x)

  list(iter = i, loglik = loglik, param = a, runtime = runtime)
}
