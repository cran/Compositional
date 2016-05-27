################################
#### One sample hypothesis testing about the mean Hotelling's T-square test
#### Tsagris Michail 8/2012
#### mtsagris@yahoo.gr
#### References: Mardia K.V., kent J.T. & Bibby J.M. (1979)
#### Multivariate Analysis p. 126. Academic Press
################################

hotel1T2 <- function(x, M, a = 0.05, R = 999, graph = FALSE) {
  ## x is the data set
  ## a is the level of significance set by default to 0.05
  ## M is the hypothesised mean
  ## R is the number of bootstrap replicates set by default to 999
  ## if R=1 no bootstrap will be implemented
  ## Bootstrap is used for the p-value

  x <- as.matrix(x)
  M <- as.vector( as.matrix(M) )
  m <- colMeans(x)  ## sample mean vector
  s <- cov(x)  ## sample covariance matrix
  n <- nrow(x)  ## sample size
  p <- ncol(x)  ## dimensionality of the data
  test <- as.vector( (n * (n - p) ) / ( (n - 1) * p ) * mahalanobis(m, M, s) )

  ## test is the test statistic
  if (R == 1) {
    pvalue <- 1 - pf(test, p, n - p)  ## p-value of the test statistic
    crit <- qf(1 - a, p, n - p)  ## critival value of the F distribution
    info <- c(test, pvalue, crit, p, n - p)
    names(info) <- c("test", "p-value", "critical", "numer df", "denom df")
    result <- list(m = m, info = info)
  }

  if (R > 1) {
    runtime <- proc.time()
    ## bootstrap calibration
    tb <- numeric(R)
    y <- x - rep( m, rep(n, p) ) + rep( M, rep(n, p) )  ## brings the data
    ## under the null hypothesis, i.e. mean vector equal to M
    for (i in 1:R) {
      b <- sample(1:n, n, replace = TRUE)
      sb <- cov(y[b, ])
      mb <- colMeans(y[b, ])
      tb[i] <- mahalanobis(mb, M, sb)
    }

    tb <- ( n * (n - p) ) / ( (n - 1) * p ) * tb
    pvalue <- ( sum(tb > test) + 1 )/(R + 1)  ## bootstrap p-value

    if ( graph == TRUE ) {
      hist(tb, xlab = "Bootstrapped test statistic", main = " ")
      abline(v = test, lty = 2, lwd = 2)  ## The dotted vertical line
      ## is the test statistic value
    }
    runtime <- proc.time() - runtime

    result <- list(m = m, pvalue = pvalue, runtime = runtime)
  }

  result
}
