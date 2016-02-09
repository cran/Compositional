################################
#### Dirichlet distribution parameters
#### Tsagris Michail 3/2012  
#### mtsagris@yahoo.gr
################################

diri.est <- function(x, type = 'mle') {
  ## x is the compositional data
  x <- as.matrix(x)  ## makes sure x is a matrix
  x <- x/rowSums(x)  ## makes sure x is compositional data
  n <- nrow(x) ## the sample size 
  ## type indicates how to estimate parameters
  ## type = 'mle' means the classical mle case
  ## type = 'prec' means to use the precision parameter phi
  ## type = 'ent' means to use the entropy for the estimation
  ## loglik is for the 'mle' type
  loglik <- function(param, x = x) {
    param = exp(param)
    f <-  -( n * lgamma( sum(param) ) - n * sum( lgamma(param) ) +
    sum( log(x) %*% (param - 1) ) ) 
    f
  }
  ## diri is for the 'prec' type
  diriphi <- function(param, x = x) {
    phi <- exp(param[1])  
    b <- c(1, exp(param[-1]) )
    b <- b / sum(b)
    f <-  -( n * lgamma(phi) - n * sum( lgamma(phi * b) ) +
    sum( log(x) %*% (phi * b - 1) ) ) 
    f  
  }
  ## entro is for the 'ent' type
  entro <- function(param) {
    f <- numeric( length(param) )
    ma <- colMeans(log(x))
    for (i in 1:length(f)) {
      f[i] <- ma[i] - digamma(param[i]) + digamma( sum(param) ) 
    }
    f
  } 
  # m <- colMeans(x)
  # down = -sum( m * colMeans( log( x %*% diag(1/m) ) ) )
  # s <- 0.5 * (ncol(x) - 1) / down  ## initial value for precision
  if (type == 'mle') {
    options(warn = -1)
    da <- nlm(loglik, colMeans(x) * 10, x = x, iterlim = 10000)
    da <- nlm(loglik, da$estimate, x = x, iterlim = 10000)
    da <- nlm(loglik, da$estimate, x = x, iterlim = 10000)
    da <- optim(da$estimate, loglik, x = x, control = list(maxit = 2000), 
    hessian = T) 
    result <- list( loglik = -da$value, param = exp( da$par ) ) 
  }
  if (type == 'prec') {
    options(warn = -1)
    da <- nlm(diriphi, c(10, colMeans(x)[-1]), x = x, iterlim = 1000)
    da <- nlm(diriphi, da$estimate, x = x, iterlim = 1000)
    da <- nlm(diriphi, da$estimate, x = x, iterlim = 1000, hessian = T)
    phi <- exp(da$estimate[1]) 
    a <- c(1, exp(da$estimate[-1]))
    a <- a/sum(a)
    result <- list( loglik = -da$minimum, phi = phi, a = a, param = phi*a ) }
  if (type == 'ent') {
    ## this requires the BB package
    da <- BB::BBsolve(runif(ncol(x), 0, 20), entro, control = 
    list(maxit = 20000, tol = 1e-10))
    da <- BB::BBsolve(da$par, entro, control = list(maxit = 20000, tol = 1e-10))
    da <- BB::BBsolve(da$par, entro, control = list(maxit = 20000, tol = 1e-10))
    param <- da$par
    lik <- n * lgamma(sum(param)) - n * sum( lgamma(param) ) + 
    sum( log(x) %*% (param - 1) )
    result <- list( loglik = lik, param = param )
  }
  result 
}