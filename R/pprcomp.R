pprcomp <- function(y, x, nterms = 3, type = "log", xnew = NULL) {

  nam <- colnames(x)
  if ( type == "alr" ) {
    x <- Compositional::alr(x)
    colnames(x) <- nam[-1]
  } else {
    x <- Rfast::Log(x)
    colnames(x) <- nam
  }

  runtime <- proc.time()
  x <- as.data.frame(x)
  nam <- colnames(x)
  p <- dim(x)[2]
  if ( is.null(nam) )  colnames(x) <- paste("X", 1:p, sep = "")

  mod <- ppr(y ~., data = x, nterms = nterms)

  if ( !is.null(xnew) ) {
    if ( type == "alr" ) {
      xnew <- Compositional::alr(xnew)
    } else  xnew <- Rfast::Log(xnew)
    xnew <- as.data.frame(xnew)
    colnames(xnew) <- nam
    est <- predict(mod, newdata = xnew)
  } else  est <- NULL

  runtime <- proc.time() - runtime
  list(runtime = runtime, mod = mod, est = est)
}
