################################
#### Multinomial or Kullback-Leibler divergence based
#### regression for compositional data
#### Tsagris Michail 2/2015
#### mtsagris@yahoo.gr
#### References: Murteira, Jose M.R. and Ramalho, Joaquim J.S. (2013)
#### Regression analysis of multivariate fractional data
#### Econometric Reviews (to appear)
################################

kl.compreg <- function(y, x, B = 1, ncores = 1, xnew = NULL) {
  ## y is dependent variable, the compositional data
  ## x is the independent variable(s)
  ## B is the number of bootstrap samples used to obtain
  ## standard errors for the betas
  ## if B==1 no bootstrap is performed and no standard errors are reported
  ## if ncores=1, then 1 processor is used, otherwise
  ## more are used (parallel computing)  n <- dim(y)[1]  ## sample size
  x <- model.matrix( ~ ., data.frame(x) )
  d <- dim(y)[2] - 1  ## dimensionality of the simplex
  n <- dim(y)[1]

   klreg <- function(para, y, x, d) {
     be <- matrix(para, byrow = TRUE, ncol = d)
     mu1 <- cbind( 1, exp(x %*% be) )
     mu <- mu1 / rowSums(mu1)
     - sum(y * log(mu), na.rm = TRUE)
   }

  ## the next lines minimize the reg function and obtain the estimated betas
  runtime <- proc.time()
  ini <- as.vector( t( coef( lm.fit(x, y[, -1]) ) ) )  ## initial values
  options (warn = -1)
  qa <- nlm(klreg, ini, y = y, x = x, d = d)
  qa <- nlm(klreg, qa$estimate, y = y, x = x, d = d)
  qa <- nlm(klreg, qa$estimate, y = y, x = x, d = d)
  be <- matrix(qa$estimate, byrow = TRUE, ncol = d)
  seb <- NULL
  runtime <- proc.time() - runtime

  if (B > 1) {
  nc <- ncores
    if (nc == 1) {
      runtime <- proc.time()
      betaboot <- matrix(nrow = B, ncol = length(ini))
      for (i in 1:B) {
        ida <- sample(1:n, n, replace = TRUE)
        yb <- y[ida, ]
        xb <- x[ida, ]
        ini <- as.vector( t( coef( lm.fit(xb, yb[, -1]) ) ) )  ## initial values
        qa <- nlm(klreg, ini, y = yb, x = xb, d = d)
        qa <- nlm(klreg, qa$estimate, y = yb, x = xb, d = d)
        qa <- nlm(klreg, qa$estimate, y = yb, x = xb, d = d)
        betaboot[i, ] <- qa$estimate
      }
      s <- Rfast::colVars(betaboot, std = TRUE)
      seb <- matrix(s, byrow = TRUE, ncol = d)
      runtime <- proc.time() - runtime

    } else {
      runtime <- proc.time()
      betaboot <- matrix(nrow = B, ncol = length(ini) )
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ww <- foreach::foreach(i = 1:B, .combine = rbind) %dopar% {
        ida <- sample(1:n, n, replace = TRUE)
          yb <- y[ida, ]
          xb <- x[ida, ]
          ini <- as.vector( t( coef( lm.fit(xb, yb[, -1]) ) ) )  ## initial values
          qa <- nlm(klreg, ini, y = yb, x = xb, d = d)
          qa <- nlm(klreg, qa$estimate, y = yb, x = xb, d = d)
          qa <- nlm(klreg, qa$estimate, y = yb, x = xb, d = d)
          betaboot[i, ] <- qa$estimate
       }
      stopCluster(cl)

      s <- Rfast::colVars(ww, std = TRUE)
      seb <- matrix(s, byrow = TRUE, ncol = d)
      runtime <- proc.time() - runtime
    }
  }

  if ( is.null(xnew) ) {
    mu <- cbind( 1, exp(x %*% be) )
    est <- mu / Rfast::rowsums(mu)
  } else {
    xnew <- model.matrix(~., data.frame(xnew) )
    mu <- cbind(1, exp(xnew %*% beta))
    est <- mu / Rfast::rowsums(mu)
  }

  rownames(be) <- colnames(x)
  if  ( !is.null(seb) ) rownames(seb) <- colnames(x)
  list(runtime = runtime, be = be, seb = seb, est = est)
}
