aeqdist.etest <- function(x, sizes, a = 1, R = 999, ms = FALSE) {

  if ( ms ) {
    x1 <- x[1:sizes[1], ]
    x2 <- x[-c(1:sizes[1]), ]

    if ( length(a) == 1 ) {
      res <- Compositional::eqdist.etest(x1, x2, R = R)
    } else {
      if ( min(x) == 0 )  a <- a[ a > 0 ]
      len <- length(a)
      res <- numeric(len)
      for ( i in 1:len )  {
        y1 <- Compositional::alfa(x1, a[i])$aff
        y2 <- Compositional::alfa(x2, a[i])$aff
        res[i] <- Compositional::eqdist.etest(y1, y2, R = R)
      }
      names(res) <- a
    }

  } else {

    if ( length(a) == 1 ) {
      x <- Compositional::alfa(x, a)$aff
      res <- energy::eqdist.etest(x, sizes, R = R)$p.value
    } else {
      if ( min(x) == 0 )  a <- a[ a > 0 ]
      len <- length(a)
      res <- numeric(len)
      for (i in 1:len) {
        y <- Compositional::alfa(x, a[i])$aff
        res[i] <- energy::eqdist.etest(y, sizes, R = R)$p.value
      }
      names(res) <- a
    }
  }

  res
}
