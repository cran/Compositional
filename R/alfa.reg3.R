alfa.reg3 <- function(y, x, a = c(-1, 1), xnew = NULL) {

  reg <- function(para, ya, ax, a, ha, d, D) {
    be <- matrix(para, ncol = d)
    zz <- cbind( 1, exp(ax %*% be) )
    ta <- rowSums(zz)
    za <- zz / ta
    ma <- ( D / a * za - 1/a ) %*% ha
    as.vector(ya - ma)
  }

  reg2 <- function(a, y, ha, x, d, D) {
    ya <- Compositional::alfa(y, a)$aff
    ax <- a * x
    ini <- as.vector( solve(crossprod(x), crossprod(x, ya) ) )
    suppressWarnings({
      mod <- minpack.lm::nls.lm( par = ini, fn = reg, ya = ya, ax = ax, a = a, ha = ha, d = d, D = D,
                                 control = minpack.lm::nls.lm.control(maxiter = 5000) )
    })
    be <- matrix(mod$par, ncol = d)
    est <- cbind( 1, exp(x %*% be) )
    est <- est/Rfast::rowsums(est)
    sum(y * log(y / est), na.rm = TRUE)
  }

  x <- model.matrix(y ~., data.frame(x) )
  D <- dim(y)[2]
  d <- D - 1  ## dimensionality of the simplex
  ha <- t( Compositional::helm(D) )
  if ( min(y) == 0 )  a <- c(0.001, 1)
  ha <- t( Compositional::helm(D) )

  runtime <- proc.time()
  a <- optimize(reg2, a, y = y, ha = ha, x = x, d = d, D = D)$minimum
  runtime <- proc.time() - runtime
  res <- Compositional::alfa.reg(y = y, x = x[, -1], a = a, xnew = xnew)
  res$runtime <- res$runtime + runtime
  res$alfa <- a
  res
}





# alfa.reg3 <- function(y, x, a = c(-1, 1), xnew = NULL, seb = FALSE) {
#
#   reg <- function(para, ya, x, a, ha, d, D) {
#     be <- matrix(para, ncol = d)
#     mu1 <- cbind( 1, exp(x %*% be) )
#     zz <- mu1^a
#     ta <- rowSums(zz)
#     za <- zz / ta
#     ma <- ( D / a * za - 1/a ) %*% ha
#     - 2 * sum( diag( crossprod(ya, ma) ) ) + sum( diag( crossprod(ma) ) )
#   }
#
#   reg2 <- function(a, y, ha, x, d, D) {
#     ya <- Compositional::alfa(y, a)$aff
#     ini <- as.vector( solve(crossprod(x), crossprod(x, ya) ) )
#     qa1 <- nlminb( ini, reg, ya = ya, x = x, a = a, ha = ha, d = d, D = D, control = list(iter.max = 5000) )
#     qa1 <- optim( qa1$par, reg, ya = ya, x = x, a = a, ha = ha, d = d, D = D, control = list(maxit = 5000) )
#     qa2 <- optim( qa1$par, reg, ya = ya, x = x, a = a, ha = ha, d = d, D = D, control = list(maxit = 5000) )
#     while (qa1$value - qa2$value > 1e-04) {
#       qa1 <- qa2
#       qa2 <- optim( qa1$par, reg, ya = ya, x = x, a = a, ha = ha, d = d, D = D, control = list(maxit = 5000), hessian = TRUE )
#     }
#     be <- matrix(qa2$par, ncol = d)
#     est <- cbind( 1, exp(x %*% be) )
#     est <- est/Rfast::rowsums(est)
#     mean(y * log(y / est), na.rm = TRUE)
#   }
#
#   x <- model.matrix(y ~., data.frame(x) )
#   D <- dim(y)[2]
#   d <- D - 1  ## dimensionality of the simplex
#   ha <- t( helm(D) )
#   if ( min(y) == 0 )  a <- c(0.001, 1)
#
#   runtime <- proc.time()
#   a <- optimize(reg2, a, y = y, ha = ha, x = x, d = d, D = D)$minimum
#   runtime <- proc.time() - runtime
#   res <- Compositional::alfa.reg(y = y, x = x[, -1], a = a, xnew = xnew, seb = seb)
#   res$runtime <- res$runtime + runtime
#   res$alfa <- a
#   res
# }
