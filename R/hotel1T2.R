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
  m <- Rfast::colmeans(x)  ## sample mean vector
  s <- Rfast::cova(x)  ## sample covariance matrix
  n <- dim(x)[1]  ## sample size
  p <- dim(x)[2]  ## dimensionality of the data
  dm <- m - M
  test <- as.vector( n * (n - p) / (n - 1) / p * dm %*% solve(s, dm) )
  ## test is the test statistic
  if (R <= 1) {
    pvalue <- pf(test, p, n - p, lower.tail = FALSE)  ## p-value of the test statistic
    crit <- qf(1 - a, p, n - p)  ## critival value of the F distribution
    info <- c(test, pvalue, crit, p, n - p)
    names(info) <- c("test", "p-value", "critical", "numer df", "denom df")
    result <- list(m = m, info = info)
  }

  if (R > 1) {
    ## bootstrap calibration
    tb <- numeric(R)
    mm <-  - m + M
    y <- Rfast::eachrow(x, mm, oper = "+") ## brings the data
    ## under the null hypothesis, i.e. mean vector equal to M
    for (i in 1:R) {
      b <- Rfast2::Sample.int(n, n, replace = TRUE)
      yb <- y[b, ]
      sb <- Rfast::cova(yb)
      mb <- Rfast::colmeans(yb)
      dmb <- mb - M
      tb[i] <- dmb %*% solve(sb, dmb)
    }
    tb <- n * (n - p) / (n - 1) / p * tb
    pvalue <- ( sum(tb > test) + 1 )/(R + 1)  ## bootstrap p-value

    if ( graph ) {
      hist(tb, xlab = "Bootstrapped test statistic", main = " ")
      abline(v = test, lty = 2, lwd = 2)  ## The dotted vertical line
      ## is the test statistic value
    }
    result <- list(m = m, pvalue = pvalue)
  }

  result
}
