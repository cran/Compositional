comp.nb <- function(xnew = NULL, x, ina, type = "beta"){

  if ( type == "beta" ) {
    mod <- Rfast2::beta.nb(xnew, x, ina)
  } else if ( type == "logitnorm" ) {
    mod <- Rfast2::logitnorm.nb(xnew, x, ina)
  } else if (type == "cauchy" ) {
    y <- Compositional::alr(x)
    ynew <- NULL
    if ( is.null(xnew) )  ynew <- Compositional::alr(xnew)
    mod <- Rfast2::cauchy.nb(xnew, x, ina)
  } else if (type == "laplace" ) {
    y <- Compositional::alr(x)
    ynew <- NULL
    if ( is.null(xnew) )  ynew <- Compositional::alr(xnew)
    mod <- Rfast2::laplace.nb(xnew, x, ina)
  } else if (type == "gamma" ) {
    y <-  - Rfast::Log(x)
    ynew <- NULL
    if ( is.null(xnew) )  ynew <-  -Rfast::Log(xnew)
    mod <- Rfast::gammanb(xnew, x, ina)
  } else if (type == "normlog" ) {
    y <-  - Rfast::Log(x)
    ynew <- NULL
    if ( is.null(xnew) )  ynew <-  - Rfast::Log(xnew)
    mod <- Rfast2::normlog.nb(xnew, x, ina)
  } else if (type == "weibull" ) {
    y <-  - Rfast::Log(x)
    ynew <- NULL
    if ( is.null(xnew) )  ynew <-  - Rfast::Log(xnew)
    mod <- Rfast2::weibull.nb(xnew, x, ina)
  }

  mod
}
