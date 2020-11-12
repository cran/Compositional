klcompreg.boot <- function(y, x, der, der2, id, b1, n, p, d, tol = 1e-07, maxiters = 50) {
  m <- Rfast::colmeans(y)
  e <- y - rep(m, rep(n, d) )
  for (i in 1:d) {
    der[id[, i]] <- Rfast::colsums( e[, i] * x )
    for (j in i:d) {
      if (i != j) {
        der2[id[, i], id[, j]] <- der2[id[, j], id[, i]] <-  - crossprod(m[i] * m[j] * x, x)
      } else  der2[id[, i], id[, i]] <- crossprod(m[i] * (1 - m[i]) * x, x)
    }
  }
  b2 <- b1 + solve(der2, der)
  k <- 2
  res <- try(
  while ( sum( abs(b2 - b1) ) > tol  &  k < maxiters) {
    k <- k + 1
    b1 <- b2
    m1 <- exp(x %*% b1)
    m <- m1 / (Rfast::rowsums(m1) + 1)
    e <- y - m
    for (i in 1:d) {
      der[id[, i]] <- Rfast::colsums( e[, i] * x )
      for (j in i:d) {
        if (i != j) {
          der2[id[, i], id[, j]] <- der2[id[, j], id[, i]] <-   - crossprod(m[, i] * m[, j] * x, x)
        } else  der2[id[, i], id[, i]] <- crossprod(m[, i] * (1 - m[, i]) * x, x)
      }
    }
    b2 <- b1 + solve(der2, der)
  },
  silent = TRUE)
  if ( class(res) == "try-error" )   b2 <- b1
   m <- cbind(1, m1)
   m <- m / Rfast::rowsums(m)
   Y <- cbind( 1 - Rfast::rowsums(y), y )
   loglik <-  - sum( cbind(Y) * log(Y/m), na.rm = TRUE )
  list(loglik = loglik, be = b2)
}
