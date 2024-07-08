lc.rq <- function(y, x, z = NULL, tau = 0.5, xnew = NULL, znew = NULL) {

  p <- dim(x)[2]

  x <- Compositional::alr(x)

  if ( is.null(z) ) {

    mod <- quantreg::rq(y~., data = as.data.frame(x), tau = tau)
    be <- c( mod$be[1], -sum(mod$be[-1]), mod$be[-1] )

    est <- NULL
    if ( !is.null(xnew) )  est <- cbind(1, log(xnew) %*% be )
    res <- list( mod = mod, be = be, est = est )

  } else {

    z <- model.matrix( y~., data=as.data.frame(z) )[, -1, drop = FALSE]
    nama <- c( nama, colnames(z) )
    d <- dim(z)[2]
    x <- cbind(x, z)
    mod <- quantreg::rq(y~., data = as.data.frame(x), tau = tau)
    be <- c( mod$be[1], -sum(mod$be[2:p]), mod$be[-1] )

    est <- NULL
    if ( !is.null(xnew)  &  !is.null(znew) ) {
      znew <- model.matrix( ~., data.frame(znew) )[, -1]
      est <- cbind(1, log(xnew), znew) %*% be
    }

    res <- list( mod = mod, be = be, est = est )

  }

  res

}
