################################
#### Regression for compositional data based on the log-ratio transformation
#### Tsagris Michail 6/2014
#### mtsagris@yahoo.gr
#### References: John Aitchison (2003)
#### The Statistical Analysis of Compositional Data p. 158-160 Blackburn Press
################################

comp.reg <- function(y, x, type = "classical", xnew = NULL) {
  ## y is dependent variable, the compositional data
  ## x is the independent variable(s)
  ## type takes two values, either 'classical' or
  ## 'spatial' for spatial median regression.
  y <- as.matrix(y)
  y <- y/rowSums(y)  ## makes sure y is compositional data
  x <- as.matrix(x)
  ## alr transformation with the first component being the base
  z <- log( y[, -1] / y[, 1] )
  if (type == "classical") {
    mod <- multivreg(z, x, plot = F, xnew = xnew)  ## classical multivariate regression
    res <- mod$suma
    di <- nrow(res[, , 1])
    beta <- seb <- matrix(nrow = di, ncol = 2)
    for (i in 1:di) {
     beta[, i] <- res[, 1, i]
     seb[, i] <- res[, 2, i]
    }
    rownames(seb) <- rownames(beta) <- rownames(res[, , 1])
    colnames(seb) <- colnames(beta) <- colnames(mod$fitted)
    est1 <- mod$est
  }
  if (type == "spatial") {
    mod <- spatmed.reg(z, x, xnew = xnew)  ## spatial median regression
    beta <- mod$beta
    seb <- mod$sb
    est1 <- mod$est
  }
  est2 <- cbind(1, exp(est1))
  est <- est2/rowSums(est2)
  list(beta = beta, seb = seb, fitted = est)
}

