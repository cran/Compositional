unitweib.reg <- function(y, x, tau = 0.5) {
 
  x <- model.matrix(y~., data = as.data.frame(x) )
  n <- dim(x)[1]
  p <- mean(y)
  ly <-  - log(y)
  sly <- sum( ly )
 
  uwr <- function(param, y, x, tau, n, ly) {
    b <- exp(param[1])
    be <- param[-1] 
    lmi <- log( 1 + exp( x %*% (-be) ) )
    - n * log(b) - n * log(tau) + sum( log(lmi) ) -
    (b - 1) * sum( log( ly/lmi ) ) - log(tau) * sum( ( ly/lmi )^b )
  }

  qa1 <- optim( rnorm(dim(x)[2] + 1), uwr, y = y, x = x, tau = tau, n = n, ly = ly, 
                control = list(maxit = 5000) )
  qa2 <- optim( qa1$par, uwr, y = y, x = x, tau = tau, n = n, ly = ly, 
                control = list(maxit = 5000) )
  while (qa1$value - qa2$value > 1e-05) {
    qa1 <- qa2
    qa2 <- optim( qa1$par, uwr, y = y, x = x, tau = tau, n = n, ly = ly, 
                  control = list(maxit = 5000), hessian = TRUE )
  }
  qa2 <- optim( qa1$par, uwr, y = y, x = x, tau = tau, n = n, ly = ly, 
                  control = list(maxit = 5000), hessian = TRUE )

  param <- c( exp(qa2$par[1]), qa2$par[-1] )
  vb <- solve(qa2$hessian)  
  info <- cbind(param, sqrt( diag(vb) ), param^2/diag(vb) )
  info <- cbind(info, pchisq(info[, 3], 1, lower.tail = FALSE))
  rownames(info) <- c( "beta", colnames(x) )
  colnames(info) <- c("Estimate", "Std. error", "Wald", "p-value")
  list(loglik = - qa2$par + sly, info = info)

}

