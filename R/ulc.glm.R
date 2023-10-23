ulc.glm <- function(y, x, z = NULL, model = "logistic", xnew = NULL, znew = NULL) {

  p <- dim(x)[2]
  x <- log(x)

  if ( is.null(z) ) {

    if ( model == "logistic" ) {
      mod <- Rfast::glm_logistic(x, y)
    } else  mod <- Rfast::glm_poisson(x, y)

    est <- NULL
    if ( !is.null(xnew) ) {
      est <- cbind(1, log(xnew) ) %*% be
      if ( model == "logistic" ) {
        est <- 1 / ( 1 + exp(-est) )
      } else  est <- exp(est)
    }

    be <- as.vector(mod$be)
    res <- list( devi = mod$devi, be = be, est = est )

  } else {

    z <- model.matrix( y~., as.data.frame(z) )[, -1]
    d <- dim(z)[2]
    x <- cbind(log(x), z)

    if ( model == "logistic" ) {
      mod <- Rfast::glm_logistic(x, y)
    } else  mod <- Rfast::glm_poisson(x, y)

    est <- NULL
    if ( !is.null(xnew)  &  !is.null(znew) ) {
      znew <- model.matrix( ~., data.frame(znew) )[, -1]
      est <- cbind(1, log(xnew), znew) %*% be
      if ( model == "logistic" ) {
        est <- 1 / ( 1 + exp(-est) )
      } else  est <- exp(est)
    }

    be <- as.vector(mod$be)
    res <- list( devi = mod$devi, be = be, est = est )

  }

  res

}
