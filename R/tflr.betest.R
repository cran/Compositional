tflr.betest <- function(y, x, B, R = 999, ncores = 1) {

  yhat <- x %*% B
  kl <- y * log(y / yhat)
  kl[ is.infinite(kl) ] <- NA
  kl <- sum(kl, na.rm = TRUE)
  n <- dim(y)[1]
  pkl <- numeric(R)

  if (ncores <= 1) {
    for (i in 1:R) {
      id <- Rfast2::Sample.int(n, n)
      pkl[i] <- Compositional::tflr(y, x[id, ])$kl
    }

  } else {
    requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    pkl <- foreach::foreach(i = 1:R, .combine = "c",
                           .packages = c("Compositional", "Rfast", "Rfast2") ) %dopar% {
      id <- Rfast2::Sample.int(n, n)
      return( Compositional::tflr(y, x[id, ])$kl )
    }
  }

  ( sum(pkl < kl) + 1 ) / (R + 1)
}
