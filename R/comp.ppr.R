### compositional projection pursuit
comp.ppr <- function(y, x, nterms = 3, type = "alr", xnew = NULL, yb = NULL ) {
  if ( is.null(yb) )  {
    if ( type == "alr" ) {
      yb <- alr(y)
    } else   yb <- alfa(y, 0, h = TRUE)$aff
  }

  runtime <- proc.time()
  x <- as.data.frame(x)
  nam <- colnames(x)
  p <- dim(x)[2]
  if ( is.null(nam) )  colnames(x) <- paste("X", 1:p, sep = "")

  mod <- ppr(yb ~., data = x, nterms = nterms)

  if ( !is.null(xnew) ) {
    xnew <- as.data.frame(xnew)
    colnames(xnew) <- nam
    est1 <- predict(mod, newdata = xnew)
    if ( type == "alr" ) {
      est <- alrinv(est1)
    } else est <- alfainv(est1, 0, h = TRUE)
  } else  est <- NULL

  runtime <- proc.time() - runtime
  list(runtime = runtime, mod = mod, est = est)
}
