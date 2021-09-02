pivot <- function(x) {
  p <- dim(x)[2]
  fac <- sqrt(p:1)
  ly <- t( Rfast::Log(x) )
  gm <- Rfast::colCumSums( ly[p:2, ] ) 
  gm <- gm[c(p - 1):1, ] / ( c(p - 1):1 ) 
  t( fac[-1] / fac[-p] * ( ly[-p, ] - gm ) )
}


## Initial function
## pivot <- function(x) {
##  p <- dim(x)[2]
##  fac <- sqrt(p:1)
##  ly <- Rfast::Log(x)
##  gm <- Rfast::colCumSums( t(ly[, p:2]) ) 
##  gm <- t( gm[c(p - 1):1, ] / ( c(p - 1):1 ) )
##  z <- ly[, -1]
##  for ( i in 1:c(p - 1) ) {
##    z[, i] <- fac[i + 1] / fac[i] * (ly[, i] - gm[, i])
##  }
##  z
## }

