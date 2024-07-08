ulc.rq2 <- function(y, x, z = NULL, tau = 0.5, xnew = NULL, znew = NULL) {

  if ( is.null(z) ) {

    X <- NULL
    for ( i in 1:length(x) )  X <- cbind(X, x[[ i ]] )
    x <- log(X)
    X <- NULL

    mod <- quantreg::rq(y~., data = as.data.frame(x), tau = tau)

    est <- NULL
    if ( !is.null(xnew) )  est <- cbind(1, log(xnew) ) %*% be

    be <- as.vector(mod$be)
    res <- list( mod = mod, be = be, est = est )

  } else {

    z <- model.matrix( y~., as.data.frame(z) )[, -1]
    X <- NULL
    for ( i in 1:length(x) )  X <- cbind(X, x[[ i ]] )
    X <- NULL
    x <- cbind(log(x), z)

    mod <- quantreg::rq(y~., data = as.data.frame(x), tau = tau)

    est <- NULL
    if ( !is.null(xnew)  &  !is.null(znew) ) {
      znew <- model.matrix( ~., data.frame(znew) )[, -1]
      est <- cbind(1, log(xnew), znew) %*% be
    }

    be <- as.vector(mod$be)
    res <- list( mod = mod, be = be, est = est )

  }

  res

}
