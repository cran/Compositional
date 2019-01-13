pow <- function(x, a) {

  if ( is.vector(a) ) {
    z <- Rfast::eachrow(x, a, oper = "^")
  } else z <- x^a

  z / Rfast::rowsums(z)
}
