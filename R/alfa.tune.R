################################
#### Profile log-likelihood for choosing the value of alpha
#### Fast way
#### Tsagris Michail 5/2013
#### References: Tsagris, M. T., Preston, S., and Wood, A. T. A. (2011).
#### A data-based power transformation for
#### compositional data. In Proceedings of the 4rth Compositional Data Analysis Workshop, Girona, Spain.
#### mtsagris@yahoo.gr
################################
alfa.tune <- function(x, B = 1, ncores = 1) {
  ## x is the compositional data
  ## x must not contain any zeros
  n <- dim(x)[1]  ## sample size
  f <- (n - 1) / n
  D <- dim(x)[2]  ## number of components
  d <- D - 1  ## dimensionality of the simplex
  ja <- sum( log(x) )  ## part of the Jacobian of the alpha transformation
  con <-  - 0.5 * n * d * log(2 * pi * f) - 0.5 * (n - 1) * d + n * (d + 0.5) * log(D)

  pa <- function(a, x) {
    trans <- Compositional::alfa(x, a)
    z <- trans$aff  ## the alpha-transformation
    -0.5 * n * log( abs( det( cov(z) ) ) ) + (a - 1) * sum( log(x) ) - D * trans$sa
  }

  if (B == 1) {
    ell <- optimize(pa, c(-1, 1), x = x, maximum = TRUE )
    aff0 <- Compositional::alfa(x, 0)
    z0 <- aff0$aff
    sa <- aff0$sa  ## part of the Jacobian determinant as well
    lik0 <- con -  n/2 * log( abs( det( Rfast::cova(z0) ) ) ) - sum( log(x) ) - D * sa
    result <- c(ell$maximum, ell$objective + con, lik0)
    names(result) <- c("best alpha", "max log-lik", "log-lik at 0")

  } else {  ## bootstrap confidence intervals
    ell <- optimize(pa, c(-1, 1), x = x, maximum = TRUE )
    ab <- numeric(B)

    if (ncores == 1) {
      runtime <- proc.time()
      for (i in 1:B) {
        ind <- sample(1:n, n, replace = TRUE)
        ab[i] <- optimize(pa, c(-1, 1), x = x[ind, ], maximum = TRUE )$maximum
      }
      runtime <- proc.time() - runtime

    } else {
      runtime <- proc.time()
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      ab <- foreach::foreach( i = 1:B, .combine = rbind, .export = c("alfa", "helm") ) %dopar% {
        ind <- sample(1:n, n, replace = TRUE)
        return( optimize(pa, c(-1, 1), x = x[ind, ], maximum = TRUE )$maximum )
      }
      parallel::stopCluster(cl)
      runtime <- proc.time() - runtime
    }

    param <- c(ell$maximum, ell$objective + con, quantile( ab, c(0.025, 0.975) ) )
    names(param)[1:2] <- c("best alpha", "max log-lik")
    hist(ab, main = "Bootstrapped alpha values", xlab = expression( paste(alpha, " values", sep = "") ) )
    abline(v = ell$maximum, col = 3)
    abline(v = mean(ab), lty = 2, col = 4)
    message <- paste("The green is the best alpha value. The blue line is the bootstrap mean value of alpha.")
    result <- list(param = param, message = message, runtime = runtime )
  }
  result
}
