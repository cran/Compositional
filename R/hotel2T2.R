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
  p <- dim(x1)[2]  ## dimensionality of the data
  n1 <- dim(x1)[1]  ## size of the first sample
  n2 <- dim(x2)[1]  ## size of the second sample
  n <- n1 + n2  ## total sample size
  xbar1 <- Rfast::colmeans(x1)  ## sample mean vector of the first sample
  xbar2 <- Rfast::colmeans(x2)  ## sample mean vector of the second sample
  dbar <- xbar2 - xbar1  ## difference of the two mean vectors
  mesoi <- rbind(xbar1, xbar2)
  rownames(mesoi) <- c("Sample 1", "Sample 2")

  if ( is.null(colnames(x1)) ) {
    colnames(mesoi) <- colnames(mesoi) <- paste("X", 1:p, sep = "")
  } else  colnames(mesoi) <- colnames(x1)
  v <- ( (n1 - 1) * Rfast::cova(x1) + (n2 - 1) * Rfast::cova(x2) )/(n - 2)
  ## v is the pooled covariance matrix
  t2 <- as.numeric( dbar %*% solve(v, dbar) )
  test <- as.vector( n1 * n2 / n * (n - p - 1) * t2 / ( (n - 2) * p ) )  ## test statistic

  if ( R <= 1 ) {
    crit <- qf(1 - a, p, n - p - 1)  ## critical value of the F distribution
    pvalue <- pf(test, p, n - p - 1, lower.tail = FALSE)  ## p-value of the test statistic
    info <- c(test, pvalue, crit, p, n - p - 1)
    names(info) <- c("test", "p-value", "critical", "numer df", "denom df")
    result <- list(mesoi = mesoi, info = info)
  }

  if (R > 1) {
    ## bootstrap calibration
    mc <- Rfast::colmeans( rbind(x1, x2) )  ## the combined sample mean vector
    ## the next two rows bring the mean vectors of the two sample equal
    ## to the combined mean and thus equal under the null hypothesis
    mc1 <- mc - xbar1
    mc2 <- mc - xbar2
    y1 <- Rfast::eachrow(x1, mc1, oper = "+")
    y2 <- Rfast::eachrow(x2, mc2, oper = "+" )
    B <- round( sqrt(R) )
    bm1 <- bm2 <- matrix(nrow = B, ncol = p)
    vb1 <- vector("list", B)
    vb2 <- vector("list", B)
    tb <- matrix(0, B, B)
    for (i in 1:B) {
      b1 <- Rfast2::Sample.int(n1, n1, replace = TRUE)
      b2 <- Rfast2::Sample.int(n2, n2, replace = TRUE)
      yb1 <- y1[b1, ]    ;   yb2 <- y2[b2, ]
      bm1[i, ] <- Rfast::colmeans(yb1)
      bm2[i, ] <- Rfast::colmeans(yb2)  ## difference of the mean vectors
      vb1[[ i ]] <- crossprod(yb1) - tcrossprod( sqrt(n1) * bm1[i, ])
      vb2[[ i ]] <- crossprod(yb2) - tcrossprod( sqrt(n2) * bm2[i, ])
      ## vb is the pooled covariance matrix
    }
    for (i in 1:B) {
      for (j in 1:B) {
        vb <- vb1[[ i ]] + vb2[[ j ]]
        db <- bm1[i, ] - bm2[j, ]
        tb[i, j] <- db %*% solve(vb, db)
      }
    }
    pvalue <- ( sum(tb > t2) + 1 )/(B^2 + 1)

    if ( graph ) {
      hist(tb, xlab = "Bootstrapped test statistic", main = " ")
      abline(v = test, lty = 2, lwd = 2)  ## The line is the test statistic
    }
    result <- list(mesoi = mesoi, pvalue = pvalue)
  }

  result
}
