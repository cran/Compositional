eqdist.etest <- function(x, y, R = 999) {
  nx <- dim(x)[1]  ;  ny <- dim(y)[1]
  n <- nx + ny
  stat <- Rfast::edist(x, y)
  z <- rbind(x, y)
  pstat <- numeric(R)
  for ( i in 1:R ) {
    id <- Rfast2::Sample.int(n, nx)
    pstat[i] <- Rfast::edist(z[id, ], z[-id, ])
  }
  ( sum( pstat >= stat ) + 1 ) / (R + 1)
}
