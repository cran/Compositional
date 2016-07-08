el.test2 <- function(y1, y2, R = 0, ncores = 1, graph = FALSE) {
  ## y1 and y2 are the two matrices containing the multivariate (or univariate data)
  ## R is the type of calibration
  ## If R = 0, the chi-square distribution is used
  ## If R = 1, the James corrected chi-square distribution is used
  ## If R = 2, the MNV F distribution is used
  ## If R > 2, bootstrap calibration is implemented

  ## the next function is for the EL
  elpa <- function(mu) {
    t1 <- emplik::el.test(y1, mu, maxit = 1000)
    t2 <- emplik::el.test(y2, mu, maxit = 1000)
    g1 <- as.numeric( t1$"-2LLR" )
    g2 <- as.numeric( t2$"-2LLR" )
    if ( round( sum(t1$wts) ) + round( sum(t2$wts) ) == n ) {
      g <- g1 + g2
    } else  g <- 1000
    g
  }

  y1 <- as.matrix(y1)
  y2 <- as.matrix(y2)
  d <- ncol(y1)  ## number of variables
  n1 <- nrow(y1)   ;   n2 <- nrow(y2)  ## sample sizes
  n <- n1 + n2  ## total sample size
  s1 <- ( (n1 - 1) / n1 ) * cov(y1)  ;  s2 <- ( (n2 - 1) / n2 ) * cov(y2)
  m1 <- colMeans(y1)    ;    m2 <- colMeans(y2)  ## mean vectors
  ## mu is the estimate of the common mean vector
  v1 <- solve(s1)    ;   v2 <- solve(s2)
  a1 <- n1 * v1    ;    a2 <-  n2 * v2
  mu1 <- solve( a1 + a2, a1 %*% m1 + a2 %*% m2 )

  runtime <- proc.time()
  apot <- nlm( elpa, mu1 )
  test <- apot$minimum
  mu <- apot$estimate
  runtime <- proc.time() - runtime

  if ( R == 0 ) {
    pvalue <- pchisq(test, d, lower.tail = FALSE)
    result <- list( test = test, dof = d, pvalue = pvalue, mu = mu, runtime = runtime,
    note = paste("Chi-square approximation") )

  } else if ( R == 1 ) {
    delta <- james(y1, y2, R = 1)$info[3]
    stat <- as.numeric( test / delta )
    pvalue <- as.numeric( pchisq(test / delta, d, lower.tail = FALSE) )
    result <- list(test = test, modif.test = stat, dof = d, pvalue = pvalue, mu = mu,
    runtime = runtime, note = paste("James corrected chi-square approximation"))

  } else if ( R == 2 ) {
    dof <- james(y1, y2, R = 2)$info[5]
    v <- dof + d - 1
    stat <- as.numeric( ( dof / (v * d) ) * test )
    pvalue <- as.numeric( pf(stat, d, dof, lower.tail = FALSE) )
    dof <- c(d, v - d + 1)
    names(dof) <- c("numer df", "denom df")
    result <- list(test = test, modif.test = stat, dof = dof, pvalue = pvalue,
    mu = mu, runtime = runtime, note = paste("F approximation"))

    ## else bootstrap calibration is implemented
  } else if (R > 2) {
    ybar1 <- colMeans(y1)
    ybar2 <- colMeans(y2)
    z1 <- y1 - rep( ybar1 - mu1, rep(n1, d) )
    z2 <- y2 - rep( ybar2 - mu1, rep(n2, d) )

    if (ncores == 1) {
      runtime <- proc.time()
      tb <- numeric(R)
      for (i in 1:R) {
        b1 <- sample(1:n1, n1, replace = TRUE)
        b2 <- sample(1:n2, n2, replace = TRUE)
        y1 <- z1[b1, ]    ;    y2 <- z2[b2, ]
        apot <- nlm( elpa, mu1 )
        tb[i] <- apot$minimum
      }
      runtime <- proc.time() - runtime

    } else {
      runtime <- proc.time()
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl) ## make the cluster

      ww <- foreach( i = 1:R, .combine = cbind, .packages = "emplik",
          .export = "el.test" ) %dopar% {
            b1 <- sample(1:n1, n1, replace = TRUE)
            b2 <- sample(1:n2, n2, replace = TRUE)
            y1 <- z1[b1, ]    ;    y2 <- z2[b2, ]
            apot <- nlm( elpa, mu1 )
            tb[i] <- apot$minimum
      }

      stopCluster(cl) ## stop the cluster
      tb <- as.vector(ww)
      runtime <- proc.time() - runtime
    }

    pvalue <- ( sum(tb > test) + 1 ) / (R + 1)
    result <- list(test = test, pvalue = pvalue, mu = mu, runtime = runtime,
    note = paste("Bootstrap calibration") )

    if (graph == TRUE) {
      hist(tb, xlab = "Bootstrapped test statistic", main = " ")
      abline(v = test, lty = 2, lwd = 2)  ## The line is the test statistic
    }
  }

  result
}
