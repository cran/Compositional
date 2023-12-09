ols.compcomp.test <- function(y, x, B = 999) {

  stat <- Compositional::ols.compcomp(y, x)$mse
  tb <- numeric(B)
  n <- dim(y)[1]
  for ( i in 1:B ) {
    id <- Rfast2::Sample.int(n, n)
    tb[i] <- Compositional::ols.compcomp(y, x[id, ])$mse
  }
  ( sum(tb > stat) + 1 ) / (B + 1)

}

