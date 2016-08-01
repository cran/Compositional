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

  y <- as.matrix(y)
  y <- y/rowSums(y)  ## makes sure y is compositional data
  x <- as.matrix(x)
  p <- ncol(x)    ;    n <- nrow(x)

  if ( p == 1 ) {
    x <- as.vector(x)
    mx <- mean(x)
    s <- sd(x)
    x <- ( x - mx ) / s

  } else {
    mx <- as.vector( Rfast::colmeans(x) )
    s <- as.vector( Rfast::colVars(x, std = TRUE) )
    x <- ( t(x) - mx ) / s  ## standardize the xnew values
    x <- t(x)
  }

  x <- as.matrix( cbind(1, x) )
  d <- ncol(y) - 1  ## dimensionality of the simplex

  if ( !is.null(xnew) ) {
    ## if the xnew is the same as the x, the classical fitted values
    ## will be returned. Otherwise, the estimated values for the
    ## new x values will be returned.
    if ( p == 1 ) {
      xnew <- as.vector(xnew)
      xnew <- ( xnew - mx ) / s

    } else {
      xnew <- as.matrix(xnew)
      xnew <- ( t(xnew) - mx ) / s  ## standardize the xnew values
      xnew <- t(xnew)
    }

    xnew <- cbind(1, xnew)
  }

  ## internal function for the alfa-regression
  reg <- function(para){
    be <- matrix(para, byrow = TRUE, ncol = d)
    mu1 <- cbind( 1, exp(x %*% be) )
    zz <- mu1^a
    ta <- rowSums(zz)
    za <- zz / ta
    za <- ( ( d + 1 ) / a ) * za - 1/a
    ma <- za %*% ha
    esa <- ya - ma
    sa <- crossprod(esa) / (n - p)
    su <- solve(sa)
    f <- ( n/2 ) * log( det(sa) ) + 0.5 * sum( esa %*% su * esa )
    f
  }

  if ( a == 0 ) {
    ya <- alfa(y, a)$aff
    mod <- comp.reg(y, x[, -1], yb = yb)
    beta <- mod$beta
    seb <- mod$seb
    runtime <- mod$runtime

  } else {
    runtime <- proc.time()

    if ( is.null(yb) ) {
      ya <- alfa(y, a)$aff
    } else {
      ya <- yb
    }

    ha <- t( helm(d + 1) )
    m0 <- numeric(d)
    ini <- as.vector( coef( lm.fit(x, ya) ) )

    qa <- nlminb( ini, reg, control = list(iter.max = 1000) )
    qa <- optim( qa$par, reg, control = list(maxit = 5000) )
    qa <- optim( qa$par, reg, control = list(maxit = 5000) )
    qa <- optim( qa$par, reg, control = list(maxit = 5000) )
    qa <- optim( qa$par, reg, control = list(maxit = 5000), hessian = TRUE )

    beta <- matrix(qa$par, byrow = TRUE, ncol = d)
    seb <- sqrt( diag( solve( qa$hessian) ) )
    seb <- matrix(seb, byrow = TRUE, ncol = d)

    runtime <- proc.time() - runtime
  }

  if ( !is.null(xnew) ) {
    mu <- cbind( 1, exp(xnew %*% beta) )
    est <- mu/rowSums(mu)
  } else {
    mu <- cbind(1, exp(x %*% beta) )
    est <- mu/rowSums(mu)
  }

  if ( is.null( colnames(x) ) ) {
    p <- ncol(x) - 1
    rownames(beta) <- c("constant", paste("X", 1:p, sep = "") )
    if ( !is.null(seb) )  rownames(seb) <- c("constant", paste("X", 1:p, sep = "") )
  } else {
    rownames(beta)  <- c("constant", colnames(x)[-1] )
    if  ( !is.null(seb) ) rownames(seb) <- c("constant", colnames(x)[-1] )
  }

  list(runtime = runtime, beta = beta, seb = seb, est = est)
}

