################################
#### MANOVA when the covariance matrices are not equal (James test)
#### Tsagris Michail 10/2014
#### mtsagris@yahoo.gr
#### References: G.S.James (1954)
#### Tests of Linear Hypothese in Univariate and Multivariate Analysis
#### when the Ratios of the Population Variances are Unknown
#### Biometrika, Vol 41, No 1/2, 19-43
################################
maovjames <- function(x, ina, a = 0.05) {
  ## x contains all the groups together
  ## a is the significance level
  ina <- as.numeric(ina)  ## the group indicator variable
  ni <- tabulate(ina)  ## the group sample sizes
  k <- max(ina)  ## the number of groups
  p <- dim(x)[2]  ## the dimensionality
  ## the objects below will be used later
  me <- mi <- W <- matrix(nrow = k, ncol = p)
  ta <- numeric(k)
  wi <- array( dim = c(p, p, k) )
  ## the next for function calculates the
  ## mean vector and covariance matrix of each group
  for (i in 1:k) {
    zi <- x[ina == i, ]
    mi[i, ] <- Rfast::colmeans( zi )
    wi[, , i] <- ni[i] * Rfast::spdinv( Rfast::cova( zi ) )
    me[i, ] <- mi[i, ] %*% wi[, , i]
  }
  W <- t( colSums( aperm(wi) ) )
  Ws <- solve(W)
  ma <- Rfast::colsums(me)
  mesi <- Ws %*% ma  ## common mean vector
  t1 <- t2 <- numeric(k)
  Ip <- diag(p)
  for (i in 1:k) {
    ta[i] <- sum( (mi[i,] - mesi) * ( wi[, , i] %*% (mi[i, ] - mesi) ) )
    exa1 <- Ip - Ws %*% wi[, , i]
    t1[i] <- sum( diag(exa1) )
    t2[i] <- sum( exa1^2 )
  }

  test <- sum(ta)  ## the test statistic
  r <- p * (k - 1)
  A <- 1 + sum( t1^2/(ni - 1) ) / 2 / r
  B <- sum( t2 / (ni - 1) + t1^2 / 2 / (ni - 1) ) / r / (r + 2)
  x2 <- qchisq(1 - a, r)
  delta <- (A + B * x2)
  twoha <- x2 * delta  ## corrected critical value of the chi-square distribution
  pvalue <- pchisq(test/delta, r, lower.tail = FALSE)  ## p-value of the test statistic
  result <- c(test, delta, twoha, pvalue)
  names(result) <- c("test", "correction", "corr.critical", "p-value")
  result
}
