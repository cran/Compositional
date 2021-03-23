lc.reg <- function(y, x, xnew = NULL) {

 x <- log(x)
 dm <- dim(x)
 n <- dm[1]  ;  p <- dm[2]
 R <- rep(1, p)   ;   ca <- 0
 xxs <- solve( crossprod(x) )
 bols <- xxs %*% crossprod(x, y)
 com <- xxs %*% R %*% solve( R %*% xxs %*% R)
 be <- bols - as.vector( com %*% R %*% bols ) - ca 
 e <- y - x %*% be
 va <- sum(e^2) / (n - p + 1)
 covbe <- ( xxs - com %*% R %*% xxs ) * va
 nama <- colnames(x)
 if ( is.null(nama) )  nama <- paste("X", 1:p, sep = "")
 names(be) <- nama
 colnames(covbe) <- rownames(covbe) <- nama
 est <- NULL
 if ( !is.null(xnew) )  est <- log(xnew) %*% be

 list(be = be, covbe = covbe, va = va, residuals = e, est = est)
}







