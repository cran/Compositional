compknn.reg <- function(y, x, xnew, k = 2:10, type = "alr", yb = NULL) {
  if ( is.null(yb) )  {
    if ( type == "alr" ) {
      yb <- alr(y)
    } else   yb <- alfa(y, 0, h = TRUE)$aff
  }
  xnew <- as.matrix(xnew)
  nu <- dim(xnew)[1]
  di <- Rfast::dista(xnew, x, k = max(k), index = TRUE, trans = FALSE, square = FALSE)
  names <- paste("k=", k, sep = "")
  est <- sapply(names, function(x) NULL)
  est1 <- matrix( nrow = nu, ncol = dim(yb)[2] )
  klen <- length(k)
  for (j in 1:klen) {
    for (i in 1:nu)  est1[i, ] <- Rfast::colmeans( yb[di[1:k[j], i], , drop = FALSE] )
    est[[ j ]] <- alrinv(est1)
  }
  est
}
