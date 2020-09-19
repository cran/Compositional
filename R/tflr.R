tflr <- function(y, x, xnew = NULL) {

  runtime <- proc.time()
  be <- codalm::codalm(y, x)
  runtime <- proc.time() - runtime
  d <- dim(y)[2]
  p <- dim(x)[2]

  if ( is.null( colnames(y) ) )  {
    colnames(be) <- paste("Y", 1:d, sep = "")
  } else colnames(be) <- colnames(y)

  if ( is.null( rownames(y) ) )  {
    rownames(be) <- paste("X", 1:p, sep = "")
  } else rownames(be) <- colnames(x)

  yhat <- x %*% be
  loglik <-  - sum( y * log(y / yhat), na.rm = TRUE )
  est <- NULL
  if ( !is.null(xnew) )  est <- xnew %*% be

  list(runtime = runtime, loglik = loglik, be = be, est = est)
}
