################################
#### Spatial median regression
#### Tsagris Michail 10/2014
#### Biman Chakraborty (2003) On multivariate quantile regression
#### Journal of Statistical Planning and Inference
#### http://www.stat.nus.edu.sg/export/sites/dsap/research/documents/tr01_2000.pdf
#### mtsagris@yahoo.gr
################################

spatmed.reg <- function(y, x, xnew = NULL, tol = 1e-07, ses = FALSE) {

  n <- dim(y)[1]
  ## the desing matrix is created
  x <- model.matrix(y ~ ., data.frame(x) )
  p <- dim(x)[2]
  d <- dim(y)[2]

  medi <- function(be) {
    be <- matrix(be, nrow = p)
    est <- x %*% be
    sum( sqrt( rowSums( (y - est)^2 ) ) )
  }

  tic <- proc.time()

  B1 <- coef( lm.fit(x,  y) )
  est <- y - x %*% B1
  ww <- sqrt( Rfast::rowsums( est^2 ) )
  z <- x / ww
  a1 <- crossprod(z, x)
  a2 <- crossprod(z, y)
  B2 <- solve(a1, a2)
  i <- 2

  while ( sum( abs(B2 - B1) ) > tol ) {
    i <- i + 1
    B1 <- B2
    est <- y - x %*% B1
    ww <- sqrt( Rfast::rowsums( est^2 ) )
    ela <- which( ww == 0 )
    z <- x / ww
    if ( length(ela) > 0 )  z[ela, ] <- 0
    a1 <- crossprod(x, z)
    a2 <- crossprod(z, y)
    B2 <- solve(a1, a2)
  }

  be <- B2
  seb <- NULL

  if ( ses ) {
    ## we use nlm and optim to obtain the standard errors
    qa <- nlm(medi, as.vector(be), iterlim = 5000, hessian = TRUE)
    seb <- sqrt( diag( solve(qa$hessian) ) )
    seb <- matrix(seb, ncol = d)

    if ( is.null(colnames(y)) ) {
      colnames(seb) <- colnames(be) <- paste("Y", 1:d, sep = "")
    } else  colnames(seb) <- colnames(be) <- colnames(y)
  }

  if ( is.null(xnew) ) {
    est <- x %*% be
  } else {
    xnew <- model.matrix( ~ ., data.frame(xnew) )
    est <- xnew %*% be
  }

  if ( is.null(colnames(y)) ) {
    colnames(be) <- paste("Y", 1:d, sep = "")
  } else  colnames(be) <- colnames(y)

  rownames(be)  <- colnames(x)
  if  ( !is.null(seb) ) rownames(seb) <- colnames(x)

  runtime <- proc.time() - tic

  list(iter = i, runtime = runtime, be = be, seb = seb, est = est)

}
