################################
#### OLS regression for compositional data
#### Tsagris Michail 10/2014
#### mtsagris@yahoo.gr
#### References: Murteira, Jose M.R. and Ramalho, Joaquim J.S. (2013)
#### Regression analysis of multivariate fractional data
#### Econometric Reviews (to appear)
################################
ols.compreg <- function(y, x, con = TRUE, B = 1, ncores = 1, xnew = NULL) {
  ## y is dependent variable, the compositional data
  ## x is the independent variable(s)
  ## B is the number of bootstrap samples used to obtain
  ## standard errors for the bes
  ## if B==1 no bootstrap is performed and no standard errors are reported
  ## if ncores=1, then 1 processor is used, otherwise
  ## more are used (parallel computing)
  olsreg <- function(para, y, x, d) {
    be <- matrix(para, byrow = TRUE, ncol = d)
    mu1 <- cbind(1, exp(x %*% be))
    mu <- mu1 / rowSums(mu1)
    sum( (y - mu)^2 )
  }

  runtime <- proc.time()
  x <- model.matrix(y ~ ., data.frame(x) )
  if ( !con )  x <- x[, -1, drop = FALSE]
  p <- dim(x)[2]
  n <- dim(y)[1]  ## sample size
  d <- dim(y)[2] - 1  ## dimensionality of the simplex
  namx <- colnames(x)
  namy <- colnames(y)
  if ( is.null( namy ) )  {
    namy <- paste("Y", 2:(d + 1), sep = "")
  } else namy <- namy[-1]

  ## the next lines minimize the reg function and obtain the estimated betas
  ini <- as.vector( t( Compositional::kl.compreg(y, x[, -1], con = con)$be ) ) ## initial values
  suppressWarnings({
    qa <- nlm(olsreg, ini, y = y, x = x, d = d)
    qa <- nlm(olsreg, qa$estimate, y = y, x = x, d = d)
    qa <- nlm(olsreg, qa$estimate, y = y, x = x, d = d)
  })
  be <- matrix(qa$estimate, byrow = TRUE, ncol = d)
  covb <- NULL
  runtime <- proc.time() - runtime

  if (B > 1) {
    nc <- ncores
    if (nc <= 1) {
      runtime <- proc.time()
      betaboot <- matrix( nrow = B, ncol = length(ini) )
      for (i in 1:B) {
        ida <- Rfast2::Sample.int(n, n, replace = TRUE)
        yb <- y[ida, ]
        xb <- x[ida, ]
        ini <- as.vector( t( Compositional::kl.compreg(yb, xb[, -1], con = con)$be ) )  ## initial values
        qa <- nlm(olsreg, ini, y = yb, x = xb, d = d)
        qa <- nlm(olsreg, qa$estimate, y = yb, x = xb, d = d)
        qa <- nlm(olsreg, qa$estimate, y = yb, x = xb, d = d)
        betaboot[i, ] <- qa$estimate
      }  ##  end  for (i in 1:B) {
      covb <- cov(betaboot)
      runtime <- proc.time() - runtime

    } else {
      runtime <- proc.time()
      requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
      betaboot <- matrix(nrow = B, ncol = length(ini) )
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      betaboot <- foreach::foreach( i = 1:B, .combine = rbind, .packages = "Rfast2",
	              .export = c( "Sample.int", "olsreg" ) ) %dopar% {
        ida <- Rfast2::Sample.int(n, n, replace = TRUE)
        yb <- y[ida, ]
        xb <- x[ida, ]
        suppressWarnings({
          ini <- as.vector( t( Compositional::kl.compreg(yb, xb[, -1], con = con)$be ) )  ## initial values
          qa <- nlm(olsreg, ini, y = yb, x = xb, d = d)
          qa <- nlm(olsreg, qa$estimate, y = yb, x = xb, d = d)
          qa <- nlm(olsreg, qa$estimate, y = yb, x = xb, d = d)
        })
        betaboot[i, ] <- qa$estimate
      }  ##  end foreach
      parallel::stopCluster(cl)
      covb <- cov(betaboot)
      runtime <- proc.time() - runtime
    }  ## end if (nc <= 1) {

    nam <- NULL
    for (i in 1:p)  nam <- c(nam, paste(namy, ":", namx[i], sep = "") )
    colnames(covb) <- rownames(covb) <- nam
  }  ## end if (B > 1) {

  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew) )
    if ( !con )  xnew <- xnew[, -1, drop = FALSE]
    mu <- cbind( 1, exp(xnew %*% beta) )
    est <- mu / Rfast::rowsums(mu)
    colnames(est) <- colnames(y)	
  }

  colnames(be) <- namy
  rownames(be) <- namx
  list(runtime = runtime, be = be, covbe = covb, est = est)
}
