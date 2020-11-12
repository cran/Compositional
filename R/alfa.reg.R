################################
#### alfa-regression
#### Tsagris Michail 11/2015
#### mtsagris@yahoo.gr
#### References: Tsagris Michail (2013)
#### Regression analysis with compositional data containing zero values
#### Chilean Journal of Statistics, 6(2): 47-57
################################
alfa.reg <- function(y, x, a, xnew = NULL, yb = NULL, seb = FALSE) {
  ## y is the compositional data (dependent variable)
  ## x is the independent variables
  ## a is the value of alpha
  ## internal function for the alfa-regression
  reg <- function(para, ya, x, ha, d, n, D) {
    be <- matrix(para, ncol = d)
    mu1 <- cbind( 1, exp(x %*% be) )
    zz <- mu1^a
    ta <- rowSums(zz)
    za <- zz / ta
    ma <- ( D / a * za - 1/a ) %*% ha
    sa <- crossprod(ya - ma)
    det(sa, log = TRUE)
  }

  D <- dim(y)[2]
  d <- D - 1  ## dimensionality of the simplex
  dm <- dim(x)
  p <- dm[2]  ;    n <- dm[1]

  if ( is.null(yb) ) {
    ya <- Compositional::alfa(y, a)$aff
  } else  ya <- yb
  x <- model.matrix(ya ~., data.frame(x) )

  if ( a == 0 ) {
    mod <- Compositional::comp.reg(y, x[, -1], yb = yb)
    be <- mod$be
    seb <- mod$seb
    runtime <- mod$runtime

  } else {
    runtime <- proc.time()

    ha <- t( helm(D) )
    ini <- as.vector( solve(crossprod(x), crossprod(x, ya) ) )
    qa1 <- nlminb( ini, reg, ya = ya, x = x, ha = ha, d = d, n = n, D = D, control = list(iter.max = 2000) )
    qa1 <- optim( qa1$par, reg, ya = ya, x = x, ha = ha, d = d, n = n, D = D, control = list(maxit = 5000) )
    qa2 <- optim( qa1$par, reg, ya = ya, x = x, ha = ha, d = d, n = n, D = D, control = list(maxit = 5000) )
    while (qa1$value - qa2$value > 1e-05) {
      qa1 <- qa2
      qa2 <- optim( qa1$par, reg, ya = ya, x = x, ha = ha, d = d, n = n, D = D, control = list(maxit = 5000), hessian = TRUE )
    }
    be <- matrix(qa2$par, ncol = d)
    runtime <- proc.time() - runtime

    if ( seb ) {
      seb <- sqrt( diag( Rfast::spdinv(qa2$hessian) ) )
      seb <- matrix(seb, ncol = d)
    } else  seb <- NULL
  }  ## end if (a == 0)

  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew) )
    est <- cbind( 1, exp(xnew %*% be) )
	  est <- est/Rfast::rowsums(est)
  }

  if ( is.null( colnames(x) ) ) {
    p <- dim(x)[2] - 1
    rownames(be) <- c("constant", paste("X", 1:p, sep = "") )
    if ( !is.null(seb) )  rownames(seb) <- c("constant", paste("X", 1:p, sep = "") )
  } else {
    rownames(be)  <- c("constant", colnames(x)[-1] )
    if ( !is.null(seb) )  rownames(seb) <- c("constant", colnames(x)[-1] )
  }

  list(runtime = runtime, be = be, seb = seb, est = est)
}





### old reg function
##  reg <- function(para, ya, x, d, n, D) {
##    be <- matrix(para, ncol = d)
##    mu1 <- cbind( 1, exp(x %*% be) )
##    zz <- mu1^a
##    ta <- rowSums(zz)
##    za <- zz / ta
##    za <- D / a * za - 1/a
##    ma <- za %*% ha
##    esa <- ya - ma
##    sa <- crossprod(esa) / (n - p)
##    su <- solve(sa)
##    n/2 * log(det(sa)) + 0.5 * sum(esa %*% su * esa)
##  }
