kl.compreg <- function(y, x, B = 1, ncores = 1, xnew = NULL, tol = 1e-07, maxiters = 50) {
  runtime <- proc.time()
  mod <- Compositional::kl.compreg2(y, x, xnew = xnew, tol = tol, maxiters = maxiters)
  if (B == 1) {
    runtime <- proc.time() - runtime
    res <-  list(runtime = runtime, iters = mod$iters, loglik = mod$loglik, be = mod$be, seb = NULL, est = mod$est)
  } else {
    if (ncores <= 1) {
      X <- model.matrix(y~., data.frame(x) )
      p <- dim(X)[2]
      Y <- y[, -1, drop = FALSE]
      dm <- dim(Y)
      n <- dm[1]    ;   d <- dm[2]
      b1 <- mod$be
      betaboot <- matrix( nrow = B, ncol = prod( dim(b1) ) )
      id <- matrix(1:c(p * d), ncol = d)
      der <- numeric(d * p)
      der2 <- matrix(0, p * d, p * d)
      for (i in 1:B) {
        ida <- sample(1:n, n, replace = TRUE)
        yb <- Y[ida, ]
        xb <- X[ida, ]
        bb <- Compositional::klcompreg.boot(yb, xb, der, der2, id, b1, n, p, d, tol = tol, maxiters = maxiters)
        betaboot[i, ] <- as.vector(bb)
      }
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      X <- model.matrix(y~., data.frame(x) )
      p <- dim(X)[2]
      Y <- y[, -1, drop = FALSE]
      dm <- dim(Y)
      n <- dm[1]    ;   d <- dm[2]
      b1 <- mod$be
      betaboot <- matrix( nrow = B, ncol = prod( dim(b1) ) )
      id <- matrix(1:c(p * d), ncol = d)
      der <- numeric(d * p)
      der2 <- matrix(0, p * d, p * d)
      betaboot <- foreach::foreach(i = 1:B, .combine = rbind, .export = "klcompreg.boot",
                  .packages = "Rfast" ) %dopar% {
        ida <- sample(1:n, n, replace = TRUE)
        yb <- Y[ida, ]
        xb <- X[ida, ]
        bb <- Compositional::klcompreg.boot(yb, xb, der, der2, id, b1, n, p, d, tol = tol, maxiters = maxiters)
        return( as.vector(bb) )
      }
      parallel::stopCluster(cl)
    }
    s <- Rfast::colVars(betaboot, std = TRUE)
    seb <- matrix(s, byrow = TRUE, ncol = d)
    runtime <- proc.time() - runtime
    res <- list(runtime = runtime, iters = mod$iters, loglik = mod$loglik, be = mod$be, seb = seb, est = mod$est)
  }
  res
}






