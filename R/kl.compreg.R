kl.compreg <- function(y, x, con = TRUE, B = 1, ncores = 1, xnew = NULL, tol = 1e-07, maxiters = 50) {

  runtime <- proc.time()
  mod <- try( Compositional::kl.compreg2(y, x, con = con, xnew = xnew, tol = tol, maxiters = maxiters), silent = TRUE )
  if ( is.infinite(mod$loglik)  |  identical( class(mod), "try-error") )  {
    x <- model.matrix(y ~ ., data.frame(x) )
    x <- x[, -1, drop = FALSE]
    if ( !con )  {
       mod <- nnet::multinom(y ~ x - 1, trace = FALSE)
    } else  mod <- nnet::multinom(y ~ x, trace = FALSE)
    be <- t( coef(mod) )
    loglik <- mod$value
    iters <- maxiters
    est <- NULL
    if ( !is.null(xnew) ) {
      xnew <- model.matrix( ~., data.frame(xnew) )
      if ( !con )  xnew <- xnew[, -1, drop = FALSE]
      mu <- cbind( 1, exp(xnew %*% be) )
      est <- mu/Rfast::rowsums(mu)
    }

  } else {
    iters <- mod$iters
    loglik <- mod$loglik
    be <- mod$be
    est <- mod$est
  }  ##  end  if ( is.infinite(loglik)  |  identical( class(mod), "try-error") )  {

  if ( B == 1 ) {
    runtime <- proc.time() - runtime
    res <-  list(runtime = runtime, iters = iters, loglik = loglik, be = be, covb = NULL, est = est)
  } else {

    if ( ncores <= 1 ) {
      X <- model.matrix( y~., data.frame(x) )
      if ( !con )  X <- X[, -1, drop = FALSE]
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
        ida <- Rfast2::Sample.int(n, n, replace = TRUE)
        yb <- Y[ida, ]
        xb <- X[ida, ]
        bb <- Compositional::klcompreg.boot(yb, xb, der, der2, id, b1, n, p, d, tol = tol, maxiters = maxiters)
        if ( is.infinite(bb$loglik)  |  identical( class(bb), "try-error") )  {
          mod <- nnet::multinom(yb ~ xb, trace = FALSE)
          bb <- t( coef(mod) )
        } else  betaboot[i, ] <- as.vector(bb$be)
      }  ##  end  for (i in 1:B) {

    } else {
      oop <- options( warn = -1 )
      on.exit( options(oop) )
      requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      X <- model.matrix(y~., data.frame(x) )
      if ( !con )  X <- X[, -1, drop = FALSE]
      p <- dim(X)[2]
      Y <- y[, -1, drop = FALSE]
      dm <- dim(Y)
      n <- dm[1]    ;   d <- dm[2]
      b1 <- mod$be
      betaboot <- matrix( nrow = B, ncol = prod( dim(b1) ) )
      id <- matrix(1:c(p * d), ncol = d)
      der <- numeric(d * p)
      der2 <- matrix(0, p * d, p * d)
      betaboot <- foreach::foreach( i = 1:B, .combine = rbind, .packages = c("Rfast", "Rfast2", "nnet"),
	              .export = c("klcompreg.boot", "multinom", "Sample.int") ) %dopar% {
        ida <- Rfast2::Sample.int(n, n, replace = TRUE)
        yb <- Y[ida, ]
        xb <- X[ida, ]
        bb <- Compositional::klcompreg.boot(yb, xb, der, der2, id, b1, n, p, d, tol = tol, maxiters = maxiters)
        if ( is.infinite(bb$loglik)  |  identical( class(bb), "try-error") )  {
          mod <- nnet::multinom(yb ~ xb, trace = FALSE)
          bb <- t( coef(mod) )
        }
        return( as.vector(bb$be) )
      }  ##  end  foreach
      parallel::stopCluster(cl)
    }  ##  end if (ncores <= 1) {
    covb <- cov(betaboot)

    namx <- colnames(X)
    namy <- colnames(y)
    if ( is.null( namy ) )  {
      namy <- paste("Y", 2:(d + 1), sep = "")
    } else namy <- namy[-1]
    nam <- NULL
    for (i in 1:p)  nam <- c(nam, paste(namy, ":", namx[i], sep = "") )
    colnames(covb) <- rownames(covb) <- nam

    runtime <- proc.time() - runtime
    res <- list(runtime = runtime, iters = iters, loglik = loglik, be = be, covbe = covb, est = est)
  }
  res
}






