propreg <- function (y, x, varb = "quasi", tol = 1e-07, maxiters = 100) {
  Rfast::prop.reg(y = y, x = x, varb = varb, tol = tol, maxiters = maxiters)
}
