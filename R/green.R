green <- function(x, theta) {
  ## x contains the compositional data
  ## theta is the power parameter, usually between -1 and 1
  D <- dim(x)[2] ## number of components
  if ( D == 1 )   x <- t(x)
  if ( theta != 0 ) {
    z <- ( x^theta - 1 ) / theta
  } else {  
    z <- Rfast::Log(x)
  }
  z
}
