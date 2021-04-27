lc.reg <- function(y, x, xnew = NULL) {

 x <- log(x)
 dm <- dim(x)
 n <- dm[1]  ;  p <- dm[2]
 x <- cbind(1, x)
 R <- c(0, rep(1, p) )  ;   ca <- 0
 xxs <- solve( crossprod(x) )
 bols <- xxs %*% crossprod(x, y)
 com <- xxs %*% R %*% solve( R %*% xxs %*% R)
 be <- as.vector( bols - com %*% R %*% bols ) - ca
 e <- y - x %*% be
 va <- sum(e^2) / (n - p + 1)
 covbe <- ( xxs - com %*% R %*% xxs ) * va
 nama <- colnames(x)
 if ( is.null(nama) )  nama <- c( "constant", paste("X", 1:p, sep = "") )
 if ( nama[1] == "" )  nama[1] <- "constant"
 names(be) <- nama
 colnames(covbe) <- rownames(covbe) <- nama
 est <- NULL
 if ( !is.null(xnew) )  est <- cbind(1, log(xnew) ) %*% be

 list(be = be, covbe = covbe, va = va, residuals = e, est = est)
}







