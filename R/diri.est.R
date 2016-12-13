################################
#### Dirichlet distribution parameters
#### Tsagris Michail 3/2012
#### mtsagris@yahoo.gr
################################

diri.est <- function(x, type = 'mle') {
  ## x is the compositional data
  ## type indicates how to estimate parameters
  ## type = 'mle' means the classical mle case
  ## type = 'prec' means to use the precision parameter phi
  ## type = 'ent' means to use the entropy for the estimation
  n <- dim(x)[1]  ## sample size
  z <- log(x)
  ## loglik is for the 'mle' type
  loglik <- function(param, z, n) {
    param <- exp(param)
    - n * lgamma( sum(param) ) + n * sum( lgamma(param) ) -
    sum( z %*% (param - 1) )
  }
  ## diri is for the 'prec' type
  diriphi <- function(param, z, n) {
   phi <- exp(param[1])
   b <- c(1, exp(param[-1]) )
   b <- b / sum(b)
   - n * lgamma(phi) + n * sum( lgamma(phi * b) ) -
   sum( z %*% (phi * b - 1) )
  }

  if (type == 'mle') {
    runtime <- proc.time()
    options(warn = -1)
    da <- nlm(loglik, Rfast::colmeans(x) * 10, z = z, n = n, iterlim = 10000)
    da <- nlm(loglik, da$estimate, z = z, n = n, iterlim = 10000)
    da <- optim(da$estimate, loglik, z = z, n = n, control = list(maxit = 2000),
    hessian = TRUE)

    runtime <- proc.time() - runtime
    result <- list( loglik = -da$value, param = exp(da$par),
    std = sqrt( diag( solve(da$hessian) ) ), runtime = runtime  )
  }

  if (type == 'prec') {
    runtime <- proc.time()
    options(warn = -1)
    da <- nlm(diriphi, c(10, Rfast::colmeans(x)[-1]), z = z, n = n, iterlim = 2000)
    da <- nlm(diriphi, da$estimate, z = z, n = n, iterlim = 2000)
    da <- optim(da$estimate, diriphi, z = z, n = n, control = list(maxit = 3000),
    hessian = TRUE)
    phi <- exp(da$par[1])
    a <- c( 1, exp(da$par[-1]) )
    a <- a / sum(a)

    runtime <- proc.time() - runtime
    result <- list( loglik = -da$value, phi = phi, a = a,
    b = phi * a, runtime = runtime )
  }

  result
}


