################################
#### Symmetric Kullback-Leibler divergence based regression for compositional data
#### Tsagris Michail 5/2015
#### mtsagris@yahoo.gr
#### References: Michail Tsagris (2015)
#### A novel, divergence based, regression for compositional data
#### Proceedings of the 28th Panhellenic Statistics Conference
################################
symkl.compreg <- function(y, x, B = 1, ncores = 1, xnew = NULL) {
  ## y is dependent variable, the compositional data
  ## x is the independent variable(s)
  ## B is the number of bootstrap samples used to obtain
  ## standard errors for the betas
  ## if B==1 no bootstrap is performed and no standard errors are reported
  ## if ncores=1, then 1 processor is used, otherwise
  ## more are used (parallel computing)
  symkl <- function(para, y, x, d){
    be <- matrix(para, byrow = TRUE, ncol = d)
    mu1 <- cbind( 1, exp(x %*% be) )
    mu <- mu1 / rowSums(mu1)
    sum( (y - mu) * log(y / mu), na.rm = TRUE )
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

  ## the next lines minimize the symkl function and obtain the estimated betas
  runtime <- proc.time()
  ini <- rnorm( d * dim(x)[2] )
  oop <- options(warn = -1)
  on.exit( options(oop) )
  qa <- nlm(symkl, ini, y = y, x = x, d = d)
  qa <- nlm(symkl, qa$estimate, y = y, x = x, d = d)
  qa <- nlm(symkl, qa$estimate, y = y, x = x, d = d)
  be <- matrix(qa$estimate, byrow = TRUE, ncol = d)
  covb <- NULL
  runtime <- proc.time() - runtime

  if (B > 1) {
  betaboot <- matrix( nrow = B, ncol = length(ini) )
  nc <- ncores
    if (nc <= 1) {
      runtime <- proc.time()
      for (i in 1:B) {
        ini <- rnorm( d * dim(x)[2] )
        ida <- sample( 1:n, n, replace = TRUE )
        yb <- y[ida, ]
        xb <- x[ida, ]
        qa <- nlm(symkl, ini, y = yb, x = xb, d = d)
        qa <- nlm(symkl, qa$estimate, y = yb, x = xb, d = d)
        qa <- nlm(symkl, qa$estimate, y = yb, x = xb, d = d)
        betaboot[i, ] <- qa$estimate
      }
      covb <- cov(betaboot)
      runtime <- proc.time() - runtime

    } else {
      runtime <- proc.time()
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      betaboot <- foreach::foreach(i = 1:B, .combine = rbind, .export = "symkl" ) %dopar% {
        ida <- sample(1:n, n, replace = TRUE)
        yb <- y[ida, ]
        xb <- x[ida, ]
        ini <- rnorm( d * dim(x)[2] )
        qa <- nlm(symkl, ini, y = yb, x = xb, d = d)
        qa <- nlm(symkl, qa$estimate, y = yb, x = xb, d = d)
        qa <- nlm(symkl, qa$estimate, y = yb, x = xb, d = d)
        return(qa$estimate)
      }  ##  end foreach
      parallel::stopCluster(cl)
      covb <- cov(betaboot)
      runtime <- proc.time() - runtime
    }  ##  end if (nc < 1)
	
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
