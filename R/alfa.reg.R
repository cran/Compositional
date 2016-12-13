################################
#### alfa-regression
#### Tsagris Michail 11/2015
#### mtsagris@yahoo.gr
#### References: Tsagris Michail (2013)
#### Regression analysis with compositional data containing zero values
#### Chilean Journal of Statistics, 6(2): 47-57
################################

alfa.reg <- function(y, x, a, xnew = NULL, yb = NULL) {
  ## y is the compositional data (dependent variable)
  ## x is the independent variables
  ## a is the value of alpha
  ## internal function for the alfa-regression
  reg <- function(para, ya, x, d, n, D) {
    be <- matrix(para, byrow = TRUE, ncol = d)
    mu1 <- cbind( 1, exp(x %*% be) )
    zz <- mu1^a
    ta <- rowSums(zz)
    za <- zz / ta
    za <- D / a * za - 1/a
    ma <- za %*% ha
    esa <- ya - ma
    sa <- crossprod(esa) / (n - p)
    su <- solve(sa)
    n/2 * log( det(sa) ) + 0.5 * sum( esa %*% su * esa )
  }

  p <- NCOL(x)   ;    n <- NROW(x)

  if ( p == 1 ) {
    x <- as.vector(x)
    mx <- mean(x)
    s <- sd(x)
    x <- ( x - mx ) / s

  } else {
    mx <- Rfast::colmeans(x)
    s <- Rfast::colVars(x, std = TRUE)
    x <- t( ( t(x) - mx ) / s )  ## standardize the xnew values
  }

  x <- cbind(1, x)

  D <- dim(y)[2]
  d <- D - 1  ## dimensionality of the simplex

  if ( !is.null(xnew) ) {
    ## if the xnew is the same as the x, the classical fitted values
    ## will be returned. Otherwise, the estimated values for the
    ## new x values will be returned.
    if ( p == 1 ) {
      xnew <- as.vector(xnew)
      xnew <- ( xnew - mx ) / s

    } else {
      xnew <- as.matrix(xnew)
      xnew <- t( ( t(xnew) - mx ) / s )  ## standardize the xnew values
    }

    xnew <- cbind(1, xnew)
  }

  if ( a == 0 ) {
    ya <- alfa(y, a)$aff
    mod <- comp.reg(y, x[, -1], yb = yb)
    be <- mod$be
    seb <- mod$seb
    runtime <- mod$runtime

  } else {
    runtime <- proc.time()

    if ( is.null(yb) ) {
      ya <- alfa(y, a)$aff
    } else  ya <- yb

    ha <- t( helm(D) )
    m0 <- numeric(d)
    ini <- as.vector( coef( lm.fit(x, ya) ) )

    qa <- nlminb( ini, reg, ya = ya, x = x, d = d, n = n, D = D, control = list(iter.max = 1000) )
    qa <- optim( qa$par, reg, ya = ya, x = x, d = d, n = n, D = D, control = list(maxit = 5000) )
    qa <- optim( qa$par, reg, ya = ya, x = x, d = d, n = n, D = D, control = list(maxit = 5000) )
    qa <- optim( qa$par, reg, ya = ya, x = x, d = d, n = n, D = D, control = list(maxit = 5000), hessian = TRUE )

    be <- matrix(qa$par, byrow = TRUE, ncol = d)
    seb <- sqrt( diag( solve( qa$hessian) ) )
    seb <- matrix(seb, byrow = TRUE, ncol = d)

    runtime <- proc.time() - runtime
  }

  if ( !is.null(xnew) ) {
    mu <- cbind( 1, exp(xnew %*% be) )
  } else   mu <- cbind(1, exp(x %*% be) )
    est <- mu / Rfast::rowsums(mu)

  if ( is.null( colnames(x) ) ) {
    p <- dim(x)[2] - 1
    rownames(be) <- c("constant", paste("X", 1:p, sep = "") )
    if ( !is.null(seb) )  rownames(seb) <- c("constant", paste("X", 1:p, sep = "") )
  } else {
    rownames(be)  <- c("constant", colnames(x)[-1] )
    if  ( !is.null(seb) ) rownames(seb) <- c("constant", colnames(x)[-1] )
  }

  list(runtime = runtime, be = be, seb = seb, est = est)
}
