lc.rq2 <- function(y, x, z = NULL, tau = 0.5, xnew = NULL, znew = NULL) {

  if ( is.null(z) ) {

    X <- NULL   ;   p <- NULL
    for ( i in 1:length(x) ) {
      X <- cbind(X, Compositional::alr( x[[ i ]] ) )
      p <- c(p, dim(x[[ i ]])[2])
    }

    mod <- quantreg::rq(y~., data = as.data.frame(X), tau = tau)

    p <- p - 1
    p <- cumsum(p)
    be <- c( mod$be[1], -sum( mod$be[1:p[1] + 1] ), mod$be[1:p[1] + 1] )
    for ( i in 2:length(p) ) {
      a <- (p[i - 1] + 1):p[i]
      be <- c(be, -sum(mod$be[a]), mod$be[a] )
    }

    est <- NULL
    if ( !is.null(xnew) ) {
      Xnew <- NULL
      for ( i in 1:length(xnew) )  Xnew <- cbind(Xnew, log( x[[ i ]] ) )
      est <- cbind(1, Xnew) %*% be
    }

    res <- list( mod = mod, be = be, est = est )

  } else {

    z <- model.matrix( y~., as.data.frame(z) )[, -1]
    d <- dim(z)[2]   ;   X <- NULL   ;   p <- NULL
    for ( i in 1:length(x) ) {
      X <- cbind(X, Compositional::alr( x[[ i ]] ) )
      p <- c(p, dim(x[[ i ]])[2])
    }
    X <- cbind(X, z)

    mod <- quantreg::rq(y~., data = as.data.frame(X), tau = tau)

    p <- p - 1
    p <- cumsum(p)
    be <- c( mod$be[1], -sum( mod$be[1:p[1] + 1] ), mod$be[1:p[1] + 1] )
    for ( i in 2:length(p) ) {
      a <- (p[i - 1] + 1):p[i]
      be <- c(be, -sum(mod$be[a]), mod$be[a] )
    }

    est <- NULL
    if ( !is.null(xnew)  &  !is.null(znew) ) {
      Xnew <- NULL
      for ( i in 1:length(x) )  Xenw <- cbind(Xnew, log( x[[ i ]] ) )
      znew <- model.matrix( ~., data.frame(znew) )[, -1]
      est <- cbind(1, xnew, znew) %*% be
    }

    res <- list( mod = mod, be = be, est = est )

  }

  res

}
