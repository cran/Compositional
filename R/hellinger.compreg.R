hellinger.compreg <- function(y, x, con = TRUE, B = 1, ncores = 1, xnew = NULL) {

  hreg <- function(para, sqy, x, d) {
    be <- matrix(para, ncol = d)
    mu1 <- cbind(1, exp(x %*% be))
    mu <- mu1 / rowSums(mu1)
    as.vector( sqy - sqrt(mu) )
  }

  x <- model.matrix(y ~ ., data.frame(x) )
  if ( !con )  x <- x[, -1, drop = FALSE]
  p <- dim(x)[2]
  n <- dim(y)[1]  ## sample size
  d <- dim(y)[2] - 1  ## dimensionality of the simplex
  namx <- colnames(x)
  namy <- colnames(y)
  if ( is.null( namy ) )  {
    namy <- paste("Y", 2:(d + 1), sep = "")
  } else namy <- namy[-1]

  runtime <- proc.time()
  ini <- as.vector( t( Compositional::kl.compreg(y, x[, -1], con = con)$be ) ) ## initial values
  mod <- minpack.lm::nls.lm( par = ini, fn = hreg, y = y, x = x, d = d,
                             control = minpack.lm::nls.lm.control(maxiter = 5000) )
  be <- matrix(mod$par, ncol = d)
  runtime <- proc.time() - runtime
  covbe <- solve(mod$hessian)

  if (B > 1) {
    nc <- ncores
    if (nc <= 1) {
      runtime <- proc.time()
      betaboot <- matrix( nrow = B, ncol = length(ini) )
      for (i in 1:B) {
        ida <- Rfast2::Sample.int(n, n, replace = TRUE)
        yb <- y[ida, ]
        xb <- x[ida, ]
        suppressWarnings({
          ini <- as.vector( t( Compositional::kl.compreg(yb, xb[, -1], con = con)$be ) )  ## initial values
          mod <- minpack.lm::nls.lm( par = ini, fn = hreg, y = y, x = x, d = d,
                                     control = minpack.lm::nls.lm.control(maxiter = 5000) )
        })
        betaboot[i, ] <- mod$par
      }  ##  end  for (i in 1:B) {
      covbe <- cov(betaboot)
      runtime <- proc.time() - runtime

    } else {
      runtime <- proc.time()
      requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
      betaboot <- matrix(nrow = B, ncol = length(ini) )
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl)
      betaboot <- foreach::foreach( i = 1:B, .combine = rbind, .packages = "Rfast2",
               .export = c( "Sample.int", "hreg" ) ) %dopar% {
               ida <- Rfast2::Sample.int(n, n, replace = TRUE)
               yb <- y[ida, ]
               xb <- x[ida, ]
               suppressWarnings({
                 ini <- as.vector( t( Compositional::kl.compreg(yb, xb[, -1], con = con)$be ) )  ## initial values
                 mod <- minpack.lm::nls.lm( par = ini, fn = hreg, y = y, x = x, d = d,
                                            control = minpack.lm::nls.lm.control(maxiter = 5000) )
               })
               betaboot[i, ] <- mod$par
        }  ##  end foreach
        parallel::stopCluster(cl)
        covbe <- cov(betaboot)
        runtime <- proc.time() - runtime
      }  ## end if (nc <= 1) {

  }  ## end if (B > 1) {

  nam <- NULL
  for (i in 1:p)  nam <- c(nam, paste(namy, ":", namx[i], sep = "") )
  colnames(covbe) <- rownames(covbe) <- nam

  est <- NULL
  if ( !is.null(xnew) ) {
  xnew <- model.matrix(~., data.frame(xnew) )
  if ( !con )  xnew <- xnew[, -1, drop = FALSE]
    mu <- cbind( 1, exp(xnew %*% beta) )
    est <- mu / Rfast::rowsums(mu)
    colnames(est) <- colnames(y)
  }

  colnames(be) <- namy
  rownames(be) <- namx
  list(runtime = runtime, be = be, covbe = covbe, est = est)
}
