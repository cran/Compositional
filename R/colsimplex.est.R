colsimplex.est <- function (x, tol = 1e-07) {
  n <- dim(x)[1]  ;  p <- dim(x)[2]
  xx <- x * (1 - x)
  simplexfun <- function(m, xx, x) sum( (x - m)^2/xx ) / ( m^2 * (1 - m)^2 )
  res <- matrix(nrow = p, ncol = 2)
  for (i in 1:p) {
    mod <- optimize(simplexfun, c(0, 1), xx = xx[, i], x = x[, i], tol = tol)
    s <- sqrt(mod$objective/n)
    res[i, ] <- c(mod$minimum, s)
  }
  loglik <-  -0.5 * n * log(2 * pi) - 1.5 * Rfast::colsums( log(xx) ) - n * log(res[, 2]) - n/2
  res <- cbind(res, loglik)
  colnames(res) <- c("mean", "sigma", "loglik")
  res
}
