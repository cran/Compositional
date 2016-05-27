################################
#### Two sample hypothesis testing about the means when the covariance
#### matrices are not equal  (James test)
#### Tsagris Michail 8/2012
#### mtsagris@yahoo.gr
#### References: G.S. James
#### Tests of Linear Hypothese in Univariate and Multivariate Analysis
#### when the Ratios of the Population Variances are Unknown (1954)
#### Biometrika, Vol 41, No 1/2, 19-43
################################
#### Two sample hypothesis testing about the means when the covariance
#### matrices are not equal  (MNV test)
#### References: Krishnamoorthy K. and Yanping Xia
#### On Selecting Tests for Equality of Two Normal Mean Vectors (2006)
#### Multivariate Behavioral Research 41(4) 533-548
################################

james <- function(y1, y2, a = 0.05, R = 999, graph = FALSE) {
  ## y1 and y2 are the two samples
  ## a is the significance level and
  ## if R==1 the James test is performed
  ## if R==2 the Nel and van der Merwe test is performed
  ## if R>2 bootstrap calculation of the p-value is performed
  ## 999 bootstrap resamples are set by default
  ## Bootstrap is used for the p-value
  ## if graph is TRUE, the bootstrap statics are plotted

  y1 <- as.matrix(y1)
  y2 <- as.matrix(y2)
  p <- ncol(y1)  ## dimensionality of the data
  n1 <- nrow(y1)  ;  n2 <- nrow(y2)  ## sample sizes
  n <- n1 + n2  ## the total sample size
  ybar1 <- colMeans(y1)  ## sample mean vector of the first sample
  ybar2 <- colMeans(y2)  ## sample mean vector of the second sample
  dbar <- ybar2 - ybar1  ## difference of the two mean vectors
  mesoi <- rbind(ybar1, ybar2)
  rownames(mesoi) <- c("Sample 1", "Sample 2")

  if ( is.null( colnames(y1) ) ) {
    colnames(mesoi) <- paste("X", 1:p, sep = "")
  } else  colnames(mesoi) <- colnames(y1)
  A1 <- cov(y1)/n1
  A2 <- cov(y2)/n2
  V <- A1 + A2  ## covariance matrix of the difference
  m0 <- numeric(p)
  test <- as.numeric( mahalanobis(dbar, m0, V) )
  b1 <- solve(V, A1)
  b2 <- solve(V, A2)
  trb1 <- sum( diag(b1) )
  trb2 <- sum( diag(b2) )

  if (R == 1) {
    ## James test
    A <- 1 + (1/(2 * p)) * ( trb1^2/(n1 - 1) + trb2^2/(n2 - 1) )
    B <- ( 1/(p * (p + 2)) ) * ( sum(b1 * b1) / (n1 - 1) +
         sum(b2 * b2)/(n2 - 1) + 0.5 * trb1 ^ 2/ (n1 - 1) + 0.5 * trb2^2/(n2 - 1) )
    x2 <- qchisq(1 - a, p)
    delta <- (A + B * x2)
    twoha <- x2 * delta  ## corrected critical value of the chi-square
    pvalue <- pchisq(test/delta, p, lower.tail = FALSE)  ## p-value of the test statistic
    info <- c(test, pvalue, delta, twoha)
    names(info) <- c("test", "p-value", "correction", "corrected.critical")
    note <- paste("James test")
    result <- list(note = note, mesoi = mesoi, info = info)
  }

  if (R == 2) {
    ## MNV test
    low <- ( sum( diag( b1 %*% b1 ) ) + trb1^2 ) / n1 +
      ( sum( diag( b2 %*% b2 ) ) + trb2^2 ) / n2
    v <- (p + p^2) / low
    test <- as.numeric( ( (v - p + 1) / (v * p) ) * test )  ## test statistic
    crit <- qf(1 - a, p, v - p + 1)  ## critical value of the F distribution
    pvalue <- pf(test, p, v - p + 1, lower.tail = FALSE)  ## p-value of the test statistic
    info <- c(test, pvalue, crit, p, v - p + 1)
    names(info) <- c("test", "p-value", "critical", "numer df", "denom df")
    note <- paste("MNV variant of James test")
    result <- list(note = note, mesoi = mesoi, info = info)
  }

  if (R > 2) {
    ## bootstrap calibration
    runtime <- proc.time()
    mc <- solve( solve(A1) + solve(A2) ) %*%
      ( solve(A1, ybar1) + solve(A2, ybar2) )
    ## mc is the combined sample mean vector
    ## the next two rows bring the mean vectors of the two sample equal
    ## to the combined mean and thus equal under the null hypothesis
    x1 <- y1 - rep( ybar1, rep(n1, p) ) + rep( mc, rep(n1, p) )
    x2 <- y2 - rep( ybar2, rep(n2, p) ) + rep( mc, rep(n2, p) )
    tb <- numeric(R)

    for (i in 1:R) {
      b1 <- sample(1:n1, n1, replace = TRUE)
      b2 <- sample(1:n2, n2, replace = TRUE)
      db <- colMeans(x1[b1, ]) - colMeans(x2[b2, ])  ## difference of the two mean vectors
      A1 <- cov(x1[b1, ]) / n1
      A2 <- cov(x2[b2, ]) / n2
      V <- A1 + A2  ## covariance matrix of the difference
      tb[i] <- as.numeric( mahalanobis(db, 0, V) )
    }

    pvalue <- ( sum(tb > test) + 1 ) / (R + 1)

    if (graph == TRUE) {
      hist(tb, xlab = "Bootstrapped test statistic", main = " ")
      abline(v = test, lty = 2, lwd = 2)  ## The line is the test statistic
    }

    note <- paste("Bootstrap calibration")
    runtime <- proc.time() - runtime
    result <- list(mesoi = mesoi, pvalue = pvalue, note = note, runtime = runtime)
  }

  result
}
