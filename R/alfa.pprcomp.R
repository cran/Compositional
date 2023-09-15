alfa.pprcomp <- function(y, x, nterms = 3, a, xnew = NULL) {

  x <- Compositional::alfa(x, a, h = FALSE)$aff
  x <- as.data.frame(x)
  nam <- colnames(x)
  p <- dim(x)[2]
  if ( is.null(nam) )  colnames(x) <- paste("X", 1:p, sep = "")

  runtime <- proc.time()

  mod <- ppr(y ~., data = x, nterms = nterms)

  if ( !is.null(xnew) ) {
    xnew <- Compositional::alfa(x, a, h = FALSE)$aff
    xnew <- as.data.frame(xnew)
    colnames(xnew) <- nam
    est <- predict(mod, newdata = xnew)
  } else  est <- NULL

  runtime <- proc.time() - runtime
  list(runtime = runtime, mod = mod, est = est)
}
