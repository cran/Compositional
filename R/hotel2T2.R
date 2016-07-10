################################
#### Two sample hypothesis testing about the means Hotelling's T-square test
#### Tsagris Michail 8/2012
#### mtsagris@yahoo.gr
#### References: Everitt Brian (2005)
#### An R and S-Plus Companion to Multivariate Analysis p. 139-140. Springer
################################

hotel2T2 <- function(x1, x2, a = 0.05, R = 999, graph = FALSE) {
  ## x1 and x2 are the two multivariate samples a is the level
  ## of significance, which by default is set to to 0.05
  ## R is the number of bootstrap replicates
  ## set by default to 999
  ## if R=1 no bootstrap will be implemented
  ## Bootstrap is used for the p-value

  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  p <- ncol(x1)  ## dimensionality of the data
  n1 <- nrow(x1)  ## size of the first sample
  n2 <- nrow(x2)  ## size of the second sample
  n <- n1 + n2  ## total sample size
  xbar1 <- colMeans(x1)  ## sample mean vector of the first sample
  xbar2 <- colMeans(x2)  ## sample mean vector of the second sample
  dbar <- xbar2 - xbar1  ## difference of the two mean vectors
  mesoi <- rbind(xbar1, xbar2)
  rownames(mesoi) <- c("Sample 1", "Sample 2")

  if ( is.null(colnames(x1)) ) {
    colnames(mesoi) <- colnames(mesoi) <- paste("X", 1:p, sep = "")
  } else  colnames(mesoi) <- colnames(x1)
  v <- ( (n1 - 1) * Rfast::cova(x1) + (n2 - 1) * Rfast::cova(x2) )/(n - 2)
  ## v is the pooled covariance matrix
  t2 <- ( n1 * n2 * (dbar %*% solve(v, dbar) ) )/n
  test <- as.vector( ( (n - p - 1) * t2 )/( (n - 2) * p ) )  ## test statistic

  if (R <= 1) {
    crit <- qf(1 - a, p, n - p - 1)  ## critical value of the F distribution
    pvalue <- pf(test, p, n - p - 1, lower.tail = FALSE)  ## p-value of the test statistic
    info <- c(test, pvalue, crit, p, n - p - 1)
    names(info) <- c("test", "p-value", "critical", "numer df", "denom df")
    result <- list(mesoi = mesoi, info = info)
  }

  if (R > 1) {
    ## bootstrap calibration
    mc <- colMeans( rbind(x1, x2) )  ## the combined sample mean vector
    ## the next two rows bring the mean vectors of the two sample equal
    ## to the combined mean and thus equal under the null hypothesis
    mc1 <-  - xbar1 + mc
    mc2 <-  - xbar2 + mc
    y1 <- x1 + rep( mc1, rep(n1, p) )
    y2 <- x2 + rep( mc2, rep(n2, p) )
    tb <- numeric(R)

    for (i in 1:R) {
      b1 <- sample(1:n1, n1, replace = TRUE)
      b2 <- sample(1:n2, n2, replace = TRUE)
      yb1 <- y1[b1, ]    ;   yb2 <- y2[b2, ]
      db <- colMeans(yb1) - colMeans(yb2)  ## difference of the mean vectors
      vb <- ( (n1 - 1) * cov(yb1) + (n2 - 1) * cov(y2) ) / (n - 2)
      ## vb is the pooled covariance matrix
      tb[i] <- ( n1 * n2 * (db %*% solve(vb, db) ) ) / n
    }
    tb <- ( (n - p - 1) * tb ) / ( (n - 2) * p )
    pvalue <- ( sum(tb > test) + 1 )/(R + 1)

    if ( graph == TRUE ) {
      hist(tb, xlab = "Bootstrapped test statistic", main = " ")
      abline(v = test, lty = 2, lwd = 2)  ## The line is the test statistic
    }

    result <- list(mesoi = mesoi, pvalue = pvalue)
  }

  result
}
