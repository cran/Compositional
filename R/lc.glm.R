lc.glm <- function(y, x, z = NULL, model = "logistic", xnew = NULL, znew = NULL) {

  p <- dim(x)[2]

  x <- Compositional::alr(x)

  if ( is.null(z) ) {

    if ( model == "logistic" ) {
      mod <- Rfast::glm_logistic(x, y)
    } else  mod <- Rfast::glm_poisson(x, y)
    be <- c( mod$be[1], -sum(mod$be[-1]), mod$be[-1] )

    est <- NULL
    if ( !is.null(xnew) ) {
      est <- cbind(1, log(xnew) %*% be )
      if ( model == "logistic" ) {
        est <- 1 / ( 1 + exp(-est) )
      } else  est <- exp(est)
    }

    res <- list( devi = mod$devi, be = be, est = est )

  } else {

    z <- model.matrix( y~., data=as.data.frame(z) )[, -1, drop = FALSE]
    nama <- c( nama, colnames(z) )
    d <- dim(z)[2]
    x <- cbind(x, z)
    if ( model == "logistic" ) {
      mod <- Rfast::glm_logistic(x, y)
    } else  mod <- Rfast::glm_poisson(x, y)
    be <- c( mod$be[1], -sum(mod$be[2:p]), mod$be[-1] )

    est <- NULL
    if ( !is.null(xnew)  &  !is.null(znew) ) {
      znew <- model.matrix( ~., data.frame(znew) )[, -1]
      est <- cbind(1, log(xnew), znew) %*% be
      if ( model == "logistic" ) {
        est <- 1 / ( 1 + exp(-est) )
      } else  est <- exp(est)
    }

    res <- list( devi = mod$devi, be = be, est = est )

  }

  res

}
