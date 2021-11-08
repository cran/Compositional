perturbation <- function(x, y, oper = "+") {

  if ( is.vector(y) ) {
    a <- Rfast::eachrow(x, y, oper = oper)
  } else {
    if (oper == "+")  a <- x * y
    if (oper == "-")  a <- x / y
  }
  
  a / Rfast::rowsums(a)
}
