lc.glm <- function(y, x, z = NULL, model = "logistic", xnew = NULL, znew = NULL) {

  p <- dim(x)[2]
  logidev <- function(be) {
    est <- x %*% be
    -2 * sum(y * est) + 2 * sum( log1p( exp(est) ) )
  }

  poisdev <- function(be) {
    est <- exp( x %*% be )
    - 2 * sum(y * est)
  }

  if ( is.null(z) ) {
    beini <- Compositional::lc.reg(y, x)$be
    x <- cbind(1, log(x))

    con <- function(be){
      f <- sum(be[-1])
      list(ceq = f, c = NULL)
    }

    if ( model == "logistic" ) {
      dev <- logidev
    } else  dev <- poisdev

    runtime <- proc.time()
    f1 <- NlcOptim::solnl( X = beini, dev, con, lb = rep(-200, p + 1), ub = rep(200, p + 1) )
    f2 <- NlcOptim::solnl( f1$par, dev, con, lb = rep(-200, p + 1), ub = rep(200, p + 1) )
    while ( f1$fn - f2$fn > 1e-04 ) {
      f1 <- f2
      f1 <- NlcOptim::solnl( f2$par, dev, con, lb = rep(-200, p + 1), ub = rep(200, p + 1) )
      f2 <- NlcOptim::solnl( f1$par, dev, con, lb = rep(-200, p + 1), ub = rep(200, p + 1) )
    }
    runtime <- proc.time() - runtime

    devi <- f2$fn + 2 * sum( y * log(y) ) * (model == "poisson")
    be <- as.vector(f2$par)
    if ( is.null( colnames(x) ) ) {
      names(be) <- c( "constant", paste("X", 1:p, sep = "") )
    } else  names(be) <- colnames(x)

    est <- NULL
    if ( !is.null(xnew) ) {
      est <- cbind(1, log(xnew) %*% be )
      if ( model == "logistic" ) {
        est <- 1 / ( 1 + exp(-est) )
      } else  est <- exp(est)
    }

    res <- list( runtime = runtime, devi = devi, be = be, est = est )

  } else {

    beini <- Compositional::lc.reg(y, x, z)$be
    z <- model.matrix( y~., as.data.frame(z) )[, -1]
    d <- dim(z)[2]
    x <- cbind(1, log(x), z)

    con <- function(be){
      f <- sum(be[2:(p + 1)])
      list(ceq = f, c = NULL)
    }

    if ( model == "logistic" ) {
      dev <- logidev
    } else  dev <- poisdev

    runtime <- proc.time()
    f1 <- NlcOptim::solnl( X = beini, dev, con, lb = rep(-200, p + d + 1), ub = rep(200, p + d + 1) )
    f2 <- NlcOptim::solnl( f1$par, dev, con, lb = rep(-200, p + d + 1), ub = rep(200, p + d + 1) )
    while ( f1$fn - f2$fn > 1e-04 ) {
      f1 <- f2
      f1 <- NlcOptim::solnl( f2$par, dev, con, lb = rep(-200, p + d + 1), ub = rep(200, p + d + 1) )
      f2 <- NlcOptim::solnl( f1$par, dev, con, lb = rep(-200, p + d + 1), ub = rep(200, p + d + 1) )
    }
    runtime <- proc.time() - runtime

    devi <- f2$fn + 2 * sum( y * log(y) ) * (model == "poisson")
    be <- as.vector(f2$par)
    nama <- colnames(x)
    if ( is.null(nama) )  nama <- c( "constant", paste("X", 1:p, sep = ""),
                                     paste("Z", 1:d, sep = "") )
    if ( nama[1] == "" )  nama[1] <- "constant"
    names(be) <- nama
    est <- NULL

    if ( !is.null(xnew)  &  !is.null(znew) ) {
      znew <- model.matrix( ~., data.frame(znew) )[, -1]
      est <- cbind(1, log(xnew), znew) %*% be
      if ( model == "logistic" ) {
        est <- 1 / ( 1 + exp(-est) )
      } else  est <- exp(est)
    }

    res <- list( runtime = runtime, devi = devi, be = be, est = est )

  }

  res

}
