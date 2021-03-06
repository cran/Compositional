kl.compreg <- function(y, x, B = 1, ncores = 1, xnew = NULL, tol = 1e-07, maxiters = 50) {

  runtime <- proc.time()
  mod <- try( Compositional::kl.compreg2(y, x, xnew = xnew, tol = tol, maxiters = maxiters), silent = TRUE )
  if ( is.infinite(mod$loglik)  |  identical( class(mod), "try-error") )  {
    mod <- nnet::multinom(y ~ x, trace = FALSE)
    be <- t( coef(mod) )
    loglik <- mod$value
    iters <- maxiters
    est <- NULL
    if ( !is.null(xnew) ) {
      xnew <- model.matrix( ~., data.frame(xnew) )
      mu <- cbind( 1, exp(xnew %*% be) )
      est <- mu/Rfast::rowsums(mu)
    }

  } else {
    iters <- mod$iters
    loglik <- mod$loglik
    be <- mod$be
    est <- mod$est
  }  ##  end  if ( is.infinite(loglik)  |  identical( class(mod), "try-error") )  {

  if (B == 1) {
    runtime <- proc.time() - runtime
    res <-  list(runtime = runtime, iters = iters, loglik = loglik, be = be, covb = NULL, est = est)
  } else {

    if (ncores <= 1) {
      X <- model.matrix( y~., data.frame(x) )
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
        if ( is.infinite(bb$loglik)  |  identical( class(bb), "try-error") )  {
          mod <- nnet::multinom(yb ~ xb, trace = FALSE)
          bb <- t( coef(mod) )
        } else  betaboot[i, ] <- as.vector(bb)
      }  ##  end  for (i in 1:B) {
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
      betaboot <- foreach::foreach(i = 1:B, .combine = rbind, .export = c("klcompreg.boot", "multinom"),
                  .packages = c("Rfast", "nnet") ) %dopar% {
        ida <- sample(1:n, n, replace = TRUE)
        yb <- Y[ida, ]
        xb <- X[ida, ]
        bb <- Compositional::klcompreg.boot(yb, xb, der, der2, id, b1, n, p, d, tol = tol, maxiters = maxiters)
        if ( is.infinite(bb$loglik)  |  identical( class(bb), "try-error") )  {
          mod <- nnet::multinom(yb ~ xb, trace = FALSE)
          bb <- t( coef(mod) )
        }
        return( as.vector(bb) )
      }  ##  end  foreach
      parallel::stopCluster(cl)
    }  ##  end if (ncores <= 1) {
    covb <- cov(betaboot)

    runtime <- proc.time() - runtime
    res <- list(runtime = runtime, iters = iters, loglik = loglik, be = be, covbe = covb, est = est)
  }
  res
}






