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
  jsreg <- function(para, y, x, d){
    be <- matrix(para, byrow = TRUE, ncol = d)
    mu1 <- cbind( 1, exp(x %*% be) )
    mu <- mu1 / rowSums(mu1)
    M <- ( mu + y ) / 2
    sum( - y * log1p(mu / y) + mu * log(mu / M), na.rm = TRUE )
  }

  n <- dim(y)[1]  ## sample size
  x <- model.matrix(y ~ ., data.frame(x) )
  p <- dim(x)[2]
  d <- dim(y)[2] - 1  ## dimensionality of the simplex
  namx <- colnames(x)
  namy <- colnames(y)
  if ( is.null( namy ) )  {
    namy <- paste("Y", 2:(d + 1), sep = "")
  } else namy <- namy[-1]

  ## the next lines minimize the kl.compreg function and obtain the estimated betas
  runtime <- proc.time()
  ini <- as.vector( t( kl.compreg(y, x[, -1, drop = FALSE])$be ) )
  oop <- options(warn = -1)
  on.exit( options(oop) )
  qa <- nlm(jsreg, ini, y = y, x = x, d = d)
  qa <- nlm(jsreg, qa$estimate, y = y, x = x, d = d)
  qa <- nlm(jsreg, qa$estimate, y = y, x = x, d = d)
  be <- matrix(qa$estimate, byrow = TRUE, ncol = d)
  covb <- NULL
  runtime <- proc.time() - runtime

  if (B > 1) {
  betaboot <- matrix( nrow = B, ncol = length(ini) )
  nc <- ncores
    if (nc <= 1) {
      runtime <- proc.time()
      for (i in 1:B) {
        ida <- sample( 1:n, n, replace = TRUE )
        yb <- y[ida, ]
        xb <- x[ida, ]
        ini <- as.vector( t( kl.compreg(yb, xb[, -1, drop = FALSE])$be ) ) ## initial values
        qa <- nlm(jsreg, ini, y = yb, x = xb, d = d)
        qa <- nlm(jsreg, qa$estimate, y = yb, x = xb, d = d)
        qa <- nlm(jsreg, qa$estimate, y = yb, x = xb, d = d)
        betaboot[i, ] <- qa$estimate
      }  ##  end for (i in 1:B) {
      covb <- cov(betaboot)
      runtime <- proc.time() - runtime

    } else {
      runtime <- proc.time()
	  oop <- options(warn = -1)
      on.exit( options(oop) )
      requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      betaboot <- foreach::foreach(i = 1:B, .combine = rbind, .export=c("jsreg", "kl.compreg") ) %dopar% {
        ida <- sample(1:n, n, replace = TRUE)
        yb <- y[ida, ]
        xb <- x[ida, ]
        ini <- as.vector( t( kl.compreg(yb, xb[, -1, drop = FALSE])$be ) ) ## initial values
        qa <- nlm(jsreg, ini, y = yb, x = xb, d = d)
        qa <- nlm(jsreg, qa$estimate, y = yb, x = xb, d = d)
        qa <- nlm(jsreg, qa$estimate, y = yb, x = xb, d = d)
        return( qa$estimate )
      }  ##  end foreach
      parallel::stopCluster(cl)
      covb <- cov(betaboot)
      runtime <- proc.time() - runtime
    }  ##  end (nc <= 1) {
	
    nam <- NULL
    for (i in 1:p)  nam <- c(nam, paste(namy, ":", namx[i], sep = "") )
    colnames(covb) <- rownames(covb) <- nam
  }  ##  end if (B > 1) {

  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew) )
    mu <- cbind( 1, exp(xnew %*% be) )
    est <- mu / Rfast::rowsums(mu)
  }  else  est <- NULL

  colnames(be) <- namy
  rownames(be) <- namx
  list(runtime = runtime, be = be, covbe = covb, est = est)
}
