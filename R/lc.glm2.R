lc.glm2 <- function(y, x, z = NULL, model = "logistic", xnew = NULL, znew = NULL) {

  if ( is.null(z) ) {

    X <- NULL   ;   p <- NULL
    for ( i in 1:length(x) ) {
      X <- cbind(X, Compositional::alr( x[[ i ]] ) )
      p <- c(p, dim(x[[ i ]])[2])
    }

    if ( model == "logistic" ) {
      mod <- Rfast::glm_logistic(X, y)
    } else  mod <- Rfast::glm_poisson(X, y)

    est <- NULL
    if ( !is.null(xnew) ) {
      est <- cbind(1, log(xnew) ) %*% be
      if ( model == "logistic" ) {
        est <- 1 / ( 1 + exp(-est) )
      } else  est <- exp(est)
    }

    p <- p - 1
    p <- cumsum(p)
    be <- c( mod$be[1], -sum( mod$be[1:p[1] + 1] ), mod$be[1:p[1] + 1] )
    for ( i in 2:length(p) ) {
      a <- (p[i - 1] + 1):p[i]
      be <- c(be, -sum(mod$be[a]), mod$be[a] )
    }

    res <- list( devi = mod$devi, be = be, est = est )

  } else {

    z <- model.matrix( y~., as.data.frame(z) )[, -1]
    d <- dim(z)[2]   ;   X <- NULL   ;   p <- NULL
    for ( i in 1:length(x) ) {
      X <- cbind(X, Compositional::alr( x[[ i ]] ) )
      p <- c(p, dim(x[[ i ]])[2])
    }
    X <- cbind(X, z)

    if ( model == "logistic" ) {
      mod <- Rfast::glm_logistic(X, y)
    } else  mod <- Rfast::glm_poisson(X, y)

    est <- NULL
    if ( !is.null(xnew)  &  !is.null(znew) ) {
      znew <- model.matrix( ~., data.frame(znew) )[, -1]
      est <- cbind(1, log(xnew), znew) %*% be
      if ( model == "logistic" ) {
        est <- 1 / ( 1 + exp(-est) )
      } else  est <- exp(est)
    }

    p <- p - 1
    p <- cumsum(p)
    be <- c( mod$be[1], -sum( mod$be[1:p[1] + 1] ), mod$be[1:p[1] + 1] )
    for ( i in 2:length(p) ) {
      a <- (p[i - 1] + 1):p[i]
      be <- c(be, -sum(mod$be[a]), mod$be[a] )
    }
    res <- list( devi = mod$devi, be = be, est = est )

  }

  res

}
