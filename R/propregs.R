propregs <- function(y, x, varb = "quasi", tol = 1e-07, logged = FALSE, maxiters = 100) {
  Rfast::prop.regs(y = y, x = x, varb = varb, tol = tol, logged = logged, maxiters = maxiters)
}
