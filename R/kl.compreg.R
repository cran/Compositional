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
  ## more are used (parallel computing)

  y <- as.matrix(y)
  y <- y / Rfast::rowsums(y)  ## makes sure y is compositional data
  n <- dim(y)[1]  ## sample size
  mat <- model.matrix(y ~ ., as.data.frame(x) )
  x <- mat[1:n, ]
  d <- dim(y)[2] - 1  ## dimensionality of the simplex
  z <- list(y = y, x = x)

   klreg <- function(para, z) {
     y <- z$y
     x <- z$x
     be <- matrix(para, byrow = TRUE, ncol = d)
     mu1 <- cbind( 1, exp(x %*% be) )
     mu <- mu1 / rowSums(mu1)
     f <-  - sum(y * log(mu), na.rm = TRUE)
     f
   }

  ## the next lines minimize the reg function and obtain the estimated betas
  ini <- as.vector( t( coef( lm.fit(x, y[, -1]) ) ) )  ## initial values

  runtime <- proc.time()
  options (warn = -1)
  qa <- nlm(klreg, ini, z = z)
  qa <- nlm(klreg, qa$estimate, z = z)
  qa <- nlm(klreg, qa$estimate, z = z)
  beta <- matrix(qa$estimate, byrow = TRUE, ncol = d)
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
        zb <- list(y = yb, x = xb)
        ini <- as.vector( t( coef( lm.fit(xb, yb[, -1]) ) ) )  ## initial values
        qa <- nlm(klreg, ini, z = zb)
        qa <- nlm(klreg, qa$estimate, z = zb)
        qa <- nlm(klreg, qa$estimate, z = zb)
        betaboot[i, ] <- qa$estimate
      }
      s <- Rfast::colVars(ww, std = TRUE)
      seb <- matrix(s, byrow = TRUE, ncol = d)
      runtime <- proc.time() - runtime

    } else {
      runtime <- proc.time()
      betaboot <- matrix(nrow = B, ncol = length(ini) )
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ww <- foreach::foreach(i = 1:B, .combine = rbind) %dopar% {
        ida <- sample(1:n, n, replace = T)
          yb <- y[ida, ]
          xb <- x[ida, ]
          zb <- list(y = yb, x = xb)
          ini <- as.vector( t( coef( lm.fit(xb, yb[, -1]) ) ) )  ## initial values
          qa <- nlm(klreg, ini, z = zb)
          qa <- nlm(klreg, qa$estimate, z = zb)
          qa <- nlm(klreg, qa$estimate, z = zb)
          betaboot[i, ] <- qa$estimate
       }
      stopCluster(cl)

      s <- Rfast::colVars(ww, std = TRUE)
      seb <- matrix(s, byrow = TRUE, ncol = d)
      runtime <- proc.time() - runtime
    }
  }

    if ( is.null(xnew) ) {
    mu <- cbind( 1, exp(x %*% beta) )
    est <- mu / Rfast::rowsums(mu)
  } else {
    xnew <- model.matrix(y ~ ., as.data.frame(xnew) )
    xnew <- xnew[1:dim(xnew)[1], ]
    mu <- cbind(1, exp(xnew %*% beta))
    est <- mu / Rfast::rowsums(mu)
  }

  if ( is.null(colnames(x)) ) {
    p <- dim(x)[2] - 1
    rownames(beta) <- c("constant", paste("X", 1:p, sep = "") )
    if ( !is.null(seb) )  rownames(seb) <- c("constant", paste("X", 1:p, sep = "") )
  } else {
    rownames(beta)  <- c("constant", colnames(x)[-1] )
    if  ( !is.null(seb) ) rownames(seb) <- c("constant", colnames(x)[-1] )
  }

  list(runtime = runtime, beta = beta, seb = seb, est = est)
}
