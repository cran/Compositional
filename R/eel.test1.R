eel.test1 <- function(x, mu, tol = 1e-06, R = 1) {

  ## x is the multivariate data
  ## xa can also be univariate data
  ## mu is the hypothesized mean

    x <- as.matrix(x)
    d <- ncol(x)
    n <- nrow(x)

    eel <- function(x, mu) {

    d <- ncol(x)
    n <- nrow(x)
    mu <- as.vector(mu)
    ## next is the root finding function
    ### step 1
    lam <- numeric(d)

    f1 <- exp( as.vector( x %*% lam ) )
    f2 <- sum(f1)
    f2a <- x * f1
    f3 <- colSums( f2a )
    f4 <- f3 / f2
    f <- f4 - mu
    der <-  - tcrossprod( f4 ) + crossprod(f2a, x) / f2
    lam <- rbind(lam, lam - solve(der, f) )
    i <- 2

    difa <- sum( abs(lam[i, ] - lam[i - 1, ] ) )
    ## step 3 and beyond

    while ( difa > tol  &  !is.na( difa ) )  {
      i <- i + 1
      f1 <- exp( as.vector( x %*% lam[i - 1, ] ) )
      f2 <- sum(f1)
      f2a <- x * f1
      f3 <- colSums( f2a )
      f4 <- f3 / f2
      f <- f4 - mu
      der <-  - tcrossprod( f4 ) + crossprod(f2a, x) / f2
      lam <- rbind(lam, lam[i - 1, ] - solve(der, f) )
      difa <- sum( abs(lam[i, ] - lam[i - 1, ] ) )
    }

    p <- f1 / f2
    stat <-  - 2 * sum( log( n * p) )
    pvalue <- pchisq(stat, d, lower.tail = FALSE)

    info <- c(stat, pvalue)
    names(info) <- c("statistic", "p-value")
    list(p = p, lambda = lam[i, ], iters = i, info = info)
    }

    runtime <- proc.time()
    res <- try( eel(x, mu), silent = TRUE )
    runtime <- proc.time() - runtime
    res$runtime <- runtime

    if ( is.na( as.numeric( res$info[1] ) ) ) {
      res$info[1] <- 1e10
      res$info[2] <- 0
      res$p <- NA
    }

    if (R > 1 ) {

      durat <- proc.time()
      test <- as.numeric( res$info[1] )

      if ( test < 1e+10 ) {
        m <- colMeans(x)
        dm <-  - m + mu
        y <- x + rep(dm, rep(n, d) )
        tb <- numeric(R)
        for ( i in 1:R ) {
          b <- sample(1:n, n, replace = TRUE)
          yb <- y[b, ]
          tb[i] <- try( eel(yb, mu), silent = TRUE )$info[1]
        }
        pvalue <- ( sum(tb > test ) + 1 ) / (R + 1)
      } else pvalue <- 1 / (R + 1)

      res$info[2] <- pvalue
      durat <- proc.time() - durat
      runtime <- runtime + durat
      res$runtime <- runtime

    }

  res

}




