################################
#### Univariate ridge regression
#### Plot to see how things go
#### Tsagris Michail 8/2015
#### mtsagris@yahoo.gr
################################

ridge.plot <- function(y, x, lambda = seq(0, 5, by = 0.1) ) {
  ## if y is a vector only
  ## x contains the independent, continuous only, variables
  ## lambda contains a grid of values of the ridge regularization parameter
  if ( all( y > 0 & y < 1 ) )  y <- log(y / ( 1 - y) ) ## logistic normal

  n <- length(y)  ## sample size
  p <- dim(x)[2]  ## dimensionality of x
  R <- length(lambda)
  be <- matrix(nrow = p, ncol = R)
  yy <- y - sum(y) / n  ## center the dependent variables
  xx <- Rfast::standardise(x)  ## standardize the independent variables
  sa <- svd(xx)
  d <- sa$d    ;    v <- t(sa$v)    ;     tu <- t(sa$u)
  d2 <- d^2    ;    A <- d * tu %*% yy

  for (i in 1:R)  be[, i] <-  crossprod( v / ( d2 + lambda[i] ), A )

  plot(lambda, be[1,], type = "l", col = 1, lty = 1,
       ylim = c( min(be), max(be) ), xlab = expression(paste(lambda, " values") ),
       ylab = "Beta coefficients")
  for (i in 2:p)  lines(lambda, be[i, ], col = i, lty = i)
}
