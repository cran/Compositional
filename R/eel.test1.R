eel.test1 <- function(x, mu, tol = 1e-06, R = 1) {
  ## x is the multivariate data
  ## xa can also be univariate data
  ## mu is the hypothesized mean
  dm <- dim(x)
  d <- dm[2]
  n <- dm[1]

   eel <- function(x, mu, n, d) {
    lam_old <- numeric(d)
    f1 <- numeric(n)
    f2 <- n
    f2a <- x
    f3 <- Rfast::colsums( f2a )
    f4 <- f3
    f <- f4 - mu
    der <-  - tcrossprod( f4 ) + crossprod(f2a, x) / f2
    lam_new <- lam_old - solve(der, f)
    i <- 2
    ## step 3 and beyond
    while ( sum( abs(lam_new - lam_old ) ) > tol  )  {
      i <- i + 1
      lam_old <- lam_new
      f1 <- exp( as.vector( x %*% lam_old ) )
      f2 <- sum(f1)
      f2a <- x * f1
      f3 <- Rfast::colsums( f2a )
      f4 <- f3 / f2
      f <- f4 - mu
      der <-  - tcrossprod( f4 ) + crossprod(f2a, x) / f2
      lam_new <- lam_old - solve(der, f)
    }
    p <- f1 / f2
    stat <-  - 2 * sum( log( n * p) )
    pvalue <- pchisq(stat, d, lower.tail = FALSE)
    info <- c(stat, pvalue)
    names(info) <- c("statistic", "p-value")
    list(p = p, lambda = lam_new, iters = i, info = info)
  }

  runtime <- proc.time()
  res <- try( eel(x, mu, n, d), silent = TRUE )
  runtime <- proc.time() - runtime
  res$runtime <- runtime

  if ( class(res) == "try-error" ) {
    res$info[1] <- 1e10
    res$info[2] <- 0
    res$p <- NA
  }

  if (R > 1 ) {
    durat <- proc.time()
    test <- as.numeric( res$info[1] )
    if ( test < 1e+10 ) {
      m <- Rfast::colmeans(x)
      dm <-  - m + mu
      y <- x + rep(dm, rep(n, d) )
      tb <- numeric(R)
      for ( i in 1:R ) {
        b <- sample(1:n, n, replace = TRUE)
        yb <- y[b, ]
        tb[i] <- try( eel(yb, mu, n, d), silent = TRUE )$info[1]
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




