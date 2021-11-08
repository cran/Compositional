alef <- function(x, a) {

  x <- as.matrix(x)
  D <- dim(x)[2]  
  if ( D == 1 ) x = t(x)
  D <- dim(x)[2] ## number of components

  if ( abs(a) < 1e-9 )  {  ## if alpha is almost zero make it zero
    xa <- log(x)
    aff <- xa - Rfast::rowmeans(xa) 
    res <- list(aff = aff) 
  } else {  
    sk <- Rfast::rowsums(x^a)
    aff <- D / a * x^a / sk - 1/a  
    res <- list(sk = sk, aff = aff) 
  }

  res
}