################################
#### Symmetric Dirichlet distribution
#### Tsagris Michail 11/2013
#### mtsagris@yahoo.gr
################################

sym.test <- function(x) {
  ## x contains the compositional data

  n <- dim(x)[1]  ## the sample size
  D <- dim(x)[2]  ## the dimensionality of the data
  szx <- sum( log(x) )

  sym <- function(a)  n * lgamma(D * a) - n * D * lgamma(a) + szx * (a - 1)

  t0 <- optimize(sym, c(0, 1000), maximum = TRUE)
  t1 <- diri.nr(x)
  a0 <- t0$maximum
  a1 <- t1$param
  h1 <- t1$loglik
  h0 <- as.numeric(t0$objective)
  test <- 2 * (h1 - h0)
  pvalue <- pchisq(test, D - 1, lower.tail = FALSE)

  if ( is.null(colnames(x)) ) {
    names(a1) <- paste("X", 1:D, sep = "")
  } else  names(a1) <- colnames(x)
  res <- c(h1, h0, test, pvalue, D - 1)

  names(res) <- c('loglik1', 'loglik0', 'test', 'pvalue', 'df')
  list(est.par = a1, one.par = a0, res = res )

}
