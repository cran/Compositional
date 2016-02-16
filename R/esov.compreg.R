################################
#### ESOV regression for compositional data
#### Tsagris Michail 5/2015
#### mtsagris@yahoo.gr
#### References: Michail Tsagris (2015)
#### A novel, divergence based, regression for compositional data
#### Proceedings of the 28th Panhellenic Statistics Conference
################################

esov.compreg <- function(y, x, B = 1000, ncores = 1, xnew = NULL) {
  ## y is dependent variable, the compositional data
  ## x is the independent variable(s)
  ## B is the number of bootstrap samples used to obtain
  ## standard errors for the betas
  ## if B==1 no bootstrap is performed and no standard errors are reported
  ## if ncores=1, then 1 processor is used, otherwise
  ## more are used (parallel computing)
  y <- as.matrix(y)
  y <- y/rowSums(y)  ## makes sure y is compositional data
  x <- as.matrix(cbind(1, x))
  d <- ncol(y) - 1  ## dimensionality of the simplex
  n <- nrow(y)  ## sample size
  z <- list(y = y, x = x)
  esovreg <- function(para, z = z){
    y <- z$y  ;  x <- z$x
    be<- matrix(para, byrow = TRUE, ncol = d)
    mu1<- cbind( 1, exp(x %*% be) )
    mu<- mu1 / rowSums(mu1)
    M <- (mu + y) / 2
    f <- sum( y * log(y / M) + mu * log(mu / M), na.rm = TRUE )
    f
  }
  ## the next lines minimize the reg function and obtain the estimated betas
  ini <- as.vector( t( coef(lm(y[, -1] ~ x[, -1])) ) )  ## initial values
  options (warn = -1)
  qa <- nlm(esovreg, ini, z = z)
  qa <- nlm(esovreg, qa$estimate, z = z)
  qa <- nlm(esovreg, qa$estimate, z = z)
  beta <- matrix(qa$estimate, byrow = TRUE, ncol = d)
  seb <- NULL
  if (B > 1) {
  betaboot <- matrix( nrow = B, ncol = length(ini) )
  nc <- ncores
    if (nc == 1) {
      for (i in 1:B) {
        ida <- sample( 1:n, n, replace = TRUE )
        yb <- y[ida, ]
        xb <- x[ida, ]
        zb <- list(y = yb, x = xb)
        ini <- as.vector( t( coef(lm(yb[, -1] ~ xb[, -1])) ) )  ## initial values
        qa <- nlm(esovreg, ini, z = zb)
        qa <- nlm(esovreg, qa$estimate, z = zb)
        qa <- nlm(esovreg, qa$estimate, z = zb)
        betaboot[i, ] <- qa$estimate
      }
     s <- apply(betaboot, 2, sd)
     seb <- matrix(s, byrow = TRUE, ncol = d)
    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ww <- foreach::foreach(i = 1:B, .combine = rbind, .export="esovreg") %dopar% {
        ida <- sample(1:n, n, replace = TRUE)
        yb <- y[ida, ]
        xb <- x[ida, ]
        zb <- list(y = yb, x = xb)
        ini <- as.vector( t( coef(lm( yb[, -1] ~ xb[, -1]) ) ) )  ## initial values
        qa <- nlm(esovreg, ini, z = zb)
        qa <- nlm(esovreg, qa$estimate, z = zb)
        qa <- nlm(esovreg, qa$estimate, z = zb)
        betaboot[i, ] <- qa$estimate
      }
      stopCluster(cl)
      s <- apply(ww, 2, sd)
      seb <- matrix(s, byrow = TRUE, ncol = d)
    }
  }
  if ( is.null(xnew) ) {
    mu <- cbind( 1, exp(x %*% beta) )
    est <- mu / rowSums(mu)
  } else {
    xnew <- cbind(1, xnew)
    xnew <- as.matrix(xnew)
    mu <- cbind(1, exp(x %*% beta))
    est <- mu / rowSums(mu)
  }
  list(beta = beta, seb = seb, est = est)
}
