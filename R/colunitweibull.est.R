colunitweibull.est <- function(x, tol = 1e-07, maxiters = 100, parallel = FALSE) {
  Rfast2::colunitweibull.mle(x, tol, maxiters, parallel)
}
