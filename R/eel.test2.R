eel.test2 <- function(y1, y2, tol = 1e-07, R = 0, graph = FALSE) {

  ## y1 and y2 are the two matrices containing the multivariate (or univariate data)
  ## R is the type of calibration
  ## If R = 0, the chi-square distribution is used
  ## If R = 1, the James corrected chi-square distribution is used
  ## If R = 2, the MNV F distribution is used
  ## If R > 2, bootstrap calibration is implemented

  x <- as.matrix(y1)
  y <- as.matrix(y2)

  eel2 <- function(x, y) {

    d <- ncol(x)
    ## next is the root finding function
    ### step 1
    lam <- numeric(d)

    fx1 <- exp( as.vector( x %*% lam ) )
    fx2 <- sum(fx1)
    fx2a <- x * fx1
    fx3 <- colSums( fx2a )
    fx4 <- fx3 / fx2

    fy1 <- exp( as.vector( - y %*% lam ) )
    fy2 <- sum(fy1)
    fy2a <- y * fy1
    fy3 <- colSums( fy2a )
    fy4 <- fy3 / fy2

    f <- fx4 - fy4
    der <-  - tcrossprod( fx4 ) + crossprod(fx2a, x) / fx2 -
            tcrossprod( fy4 ) + crossprod(fy2a, y) / fy2

    lam <- rbind( lam, lam - solve(der, f) )
    i <- 2
    difa <- sum( abs(lam[i, ] - lam[i - 1, ] ) )

    ## step 3 and beyond

    while ( difa > tol  &  !is.na( difa ) )  {
      i <- i + 1
      fx1 <- exp( as.vector( x %*% lam[i - 1, ] ) )
      fx2 <- sum(fx1)
      fx2a <- x * fx1
      fx3 <- colSums( fx2a )
      fx4 <- fx3 / fx2

      fy1 <- exp( as.vector( - y %*% lam[i - 1, ] ) )
      fy2 <- sum(fy1)
      fy2a <- y * fy1
      fy3 <- colSums( fy2a )
      fy4 <- fy3 / fy2

      f <- fx4 - fy4
      der <-  - tcrossprod( fx4 ) + crossprod(fx2a, x) / fx2 -
            tcrossprod( fy4 ) + crossprod(fy2a, y) / fy2

      lam <- rbind(lam, lam[i - 1, ] - solve(der, f) )
      difa <- sum( abs(lam[i, ] - lam[i - 1, ] ) )
    }


      p1 <- fx1 / fx2
      p2 <- fy1 / fy2
      n1 <- nrow(x)   ;   n2 <- nrow(y)
      stat <-  - 2 * sum( log( n1 * p1) ) - 2 * sum( log(n2 * p2) )
      pvalue <- pchisq(stat, d, lower.tail = FALSE)

      info <- c(stat, pvalue, d)
      names(info) <- c("statistic", "p-value", "degrees of freedom")
      list(p1 = p1, p2 = p2, lambda = lam[i, ], iters = i, info = info)
    }

  runtime <- proc.time()
  res <- try( eel2(x, y), silent = TRUE )
  runtime <- proc.time() - runtime
  res$runtime <- runtime
  res$note <- paste("Chi-square approximation")

  if ( is.na( as.numeric( res$info[1] ) ) ) {
    res$info[1] <- 1e10
    res$info[2] <- 0
    res$p1 <- NA
    res$p2 <- NA
  }

  if ( R == 0 || res$info[1] == 1e10 ) {
    res <- res

  } else if ( R == 1 ) {
    delta <- james(y1, y2, R = 1)$info[3]
    stat <- as.numeric( test / delta )
    pvalue <- as.numeric( pchisq(test / delta, d, lower.tail = FALSE) )
    result$info[1] <- stat
    result$info[2] <- pvalue
    res$note <- paste("James corrected chi-square approximation")


  } else if ( R == 2 ) {
    dof <- james(y1, y2, R = 2)$info[5]
    v <- dof + d - 1
    stat <- as.numeric( ( dof / (v * d) ) * test )
    pvalue <- as.numeric( pf(stat, d, dof, lower.tail = FALSE) )
    dof <- c(d, v - d + 1)
    result$info <- c(stat, pvalue, dof)
    names(result$info) <- c("statistic", "p-value", "numer df", "denom df")
    res$note <- paste("F approximation")

    ## else bootstrap calibration is implemented
    if ( R > 2 ) {

      durat <- proc.time()
      test <- as.numeric( res$info[1] )

      if ( test < 1e+10 ) {
        m1 <- colMeans(x)
        m2 <- colMeans(y)
        d <- ncol(x)
        n1 <- nrow(x)   ;   n2 <- nrow(y)
        s1 <- ( (n1 - 1) / n1 ) * cov(x)
        s2 <- ( (n2 - 1) / n2 ) * cov(y)
        v1 <- solve(s1)   ;  v2 <- solve(s2)
        a1 <- n1 * v1   ;   a2 <-  n2 * v2
        ## mu is the estimate of the common mean vector
        mu <- solve( a1 + a2, a1 %*% m1 + a2 %*% m2 )
        zx <- x - rep( m1, rep(n1, d) ) + rep( mu, rep(n1, d) )
        zy <- y - rep( m2, rep(n2, d) ) + rep( mu, rep(n2, d) )

        tb <- numeric(R)
        for ( i in 1:R ) {
          b1 <- sample(1:n1, n1, replace = TRUE)
          b2 <- sample(1:n2, n2, replace = TRUE)
          y1 <- zx[b1, ]   ;   y2 <- zy[b2, ]
          tb[i] <- try( eel2(y1, y2), silent = TRUE )$info[1]
        }
        pvalue <- ( sum(tb > test ) + 1 ) / (R + 1)

        if (graph == TRUE) {
          hist(tb, xlab = "Bootstrapped test statistic", main = " ")
          abline(v = test, lty = 2, lwd = 2)  ## The line is the test statistic
        }

      } else pvalue <- 1 / (R + 1)

      res$info[2] <- pvalue
      durat <- proc.time() - durat
      runtime <- runtime + durat
      res$runtime <- runtime
      res$note <- paste("Bootstrap calibration")

    }
  }

  res

}




