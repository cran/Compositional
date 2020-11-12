ols.prop.reg <- function(y, x, cov = FALSE, tol = 1e-07, maxiters = 100){
  Rfast2::propols.reg(y = y, x = x, cov = cov, tol = tol, maxiters = maxiters)
}
