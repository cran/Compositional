aeqdist.etest <- function(x, sizes, a = 1, R = 999) {
  if ( length(a) == 1 ) {
    x <- Compositional::alfa(x, a)$aff
    res <- energy::eqdist.etest(x, sizes, R = R)$p.value
  } else {
    if ( min(x) == 0 )  a = a[ a > 0 ]
    len <- length(a)
    res <- numeric(len)
    for (i in 1:len) {
      y <- Compositional::alfa(x, a[i])$aff
      res[i] <- energy::eqdist.etest(y, sizes, R = R)$p.value
    }
    names(res) <- a
  }
  res
}
