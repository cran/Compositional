################################
#### Jensen-Shannon divergence based regression for compositional data
#### Tsagris Michail 5/2015
#### mtsagris@yahoo.gr
#### References: Michail Tsagris (2015)
#### A novel, divergence based, regression for compositional data
#### Proceedings of the 28th Panhellenic Statistics Conference
################################

js.compreg <- function(y, x, B = 1, ncores = 1, xnew = NULL) {
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

  jsreg <- function(para, z = z){
    y <- z$y   ;   x <- z$x
    be <- matrix(para, byrow = TRUE, ncol = d)
    mu1 <- cbind( 1, exp(x %*% be) )
    mu <- mu1 / rowSums(mu1)
    M <- ( mu + y ) / 2
    f <- sum( - y * log(1 + mu / y) + mu * log(mu / M), na.rm = TRUE )
    f
  }

  ## the next lines minimize the kl.compreg function and obtain the estimated betas
  ini <- as.vector( t( kl.compreg(y, x[, -1])$beta ) )

  runtime <- proc.time()
  options (warn = -1)
  qa <- nlm(jsreg, ini, z = z)
  qa <- nlm(jsreg, qa$estimate, z = z)
  qa <- nlm(jsreg, qa$estimate, z = z)
  beta <- matrix(qa$estimate, byrow = TRUE, ncol = d)
  seb <- NULL
  runtime <- proc.time() - runtime

  if (B > 1) {
  betaboot <- matrix( nrow = B, ncol = length(ini) )
  nc <- ncores
    if (nc == 1) {
      runtime <- proc.time()
      for (i in 1:B) {
        ida <- sample( 1:n, n, replace = TRUE )
        yb <- y[ida, ]
        xb <- x[ida, ]
        zb <- list(y = yb, x = xb)
        ini <- as.vector( t( kl.compreg(yb, xb[, -1])$beta ) ) ## initial values
        qa <- nlm(jsreg, ini, z = zb)
        qa <- nlm(jsreg, qa$estimate, z = zb)
        qa <- nlm(jsreg, qa$estimate, z = zb)
        betaboot[i, ] <- qa$estimate
      }
      s <- Rfast::colVars(ww, std = TRUE)
      seb <- matrix(s, byrow = TRUE, ncol = d)
      runtime <- proc.time() - runtime

    } else {
      runtime <- proc.time()
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ww <- foreach::foreach(i = 1:B, .combine = rbind, .export="jsreg") %dopar% {
        ida <- sample(1:n, n, replace = TRUE)
        yb <- y[ida, ]
        xb <- x[ida, ]
        zb <- list(y = yb, x = xb)
        ini <- as.vector( t( kl.compreg(yb, xb[, -1])$beta ) ) ## initial values
        qa <- nlm(jsreg, ini, z = zb)
        qa <- nlm(jsreg, qa$estimate, z = zb)
        qa <- nlm(jsreg, qa$estimate, z = zb)
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
