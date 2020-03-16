################################
#### Regression for compositional data based on the log-ratio transformation
#### Tsagris Michail 6/2014
#### mtsagris@yahoo.gr
#### References: John Aitchison (2003)
#### The Statistical Analysis of Compositional Data p. 158-160 Blackburn Press
################################
comp.reg <- function(y, x, type = "classical", xnew = NULL, yb = NULL) {
  ## y is dependent variable, the compositional data
  ## x is the independent variable(s)
  ## type takes three values, either 'classical' or
  ## 'spatial' for spatial median regression.
  ## alr transformation with the first component being the base
  if ( is.null(yb) )  {
    z <- log( y[, -1] / y[, 1] )
  } else  z <- yb

  if (type == "lmfit") {
    runtime <- proc.time()
    x <- model.matrix(z ~ ., data.frame(x) )
    be <- solve( crossprod(x), crossprod(x, z) )
    if  ( !is.null(xnew) )  {
      xnew <- model.matrix( ~ ., data.frame(xnew) )
      est1 <- xnew %*% be
    } else  est1 <- NULL
    seb <- NULL
    runtime <- proc.time() - runtime

  } else if (type == "classical") {
    runtime <- proc.time()
    mod <- Compositional::multivreg(z, x, plot = FALSE, xnew = xnew)  ## classical multivariate regression
    res <- mod$suma
    di <- ncol(z)
    be <- seb <- matrix(nrow = NCOL(x) + 1, ncol = di)
    for (i in 1:di) {
      be[, i] <- res[, 1, i]
      seb[, i] <- res[, 2, i]
    }
    rownames(seb) <- rownames(be) <- rownames(res[, , 1])
    colnames(seb) <- colnames(be) <- colnames(mod$fitted)
    est1 <- mod$est
    runtime <- proc.time() - runtime

  } else if (type == "spatial") {
    mod <- Compositional::spatmed.reg(z, x, xnew = xnew)  ## spatial median regression
    be <- mod$be
    seb <- mod$seb
    est1 <- mod$est
    runtime <- mod$runtime
  }

  est <- NULL
  if ( !is.null(est1) ) {
    est2 <- cbind( 1, exp(est1) )
    est <- est2 / Rfast::rowsums(est2)
  }
  list(runtime = runtime, be = be, seb = seb, est = est)
}
