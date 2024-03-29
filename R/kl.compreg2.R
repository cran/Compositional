kl.compreg2 <- function(y, x, con = TRUE, xnew = NULL, tol = 1e-07, maxiters = 50) {
  X <- model.matrix( y~., data.frame(x) )
  if ( !con )  X <- X[, -1, drop = FALSE]
  p <- dim(X)[2]
  Y <- y[, -1, drop = FALSE]
  dm <- dim(Y)
  n <- dm[1]
  d <- dm[2]
  m <- Rfast::colmeans(Y)
  b0 <- Rfast::Log(m / (1 - m) )
  b1 <- matrix( c(b0, numeric(p * d - d) ), nrow = p, ncol = d, byrow = TRUE)
  e <- Y - rep(m, rep(n, d) )
  id <- matrix(1:c(p * d), ncol = d)
  der <- numeric(d * p)
  der2 <- matrix(0, p * d, p * d)
  for (i in 1:d) {
    der[id[, i]] <- Rfast::eachcol.apply(X, e[, i])
    for (j in i:d) {
      if (i != j) {
        der2[id[, i], id[, j]] <- der2[id[, j], id[, i]] <-  - crossprod(m[i] * m[j] * X, X)
      } else  der2[id[, i], id[, i]] <- crossprod(m[i] * (1 - m[i]) * X, X)
    }
  }
  b2 <- b1 + solve(der2, der)

  k <- 2
  res <- try(
  while ( sum( abs(b2 - b1) ) > tol  &  k < maxiters) {
    k <- k + 1
    b1 <- b2
    m1 <- exp(X %*% b1)
    m <- m1 / (Rfast::rowsums(m1) + 1)
    e <- Y - m
    for (i in 1:d) {
      der[id[, i]] <- Rfast::eachcol.apply(X, e[, i])
      for (j in i:d) {
        if (i != j) {
          der2[id[, i], id[, j]] <- der2[id[, j], id[, i]] <-   - crossprod(m[, i] * m[, j] * X, X)
        } else  der2[id[, i], id[, i]] <- crossprod(m[, i] * (1 - m[, i]) * X, X)
      }
    }
    b2 <- b1 + solve(der2, der)
  },
  silent = TRUE)
  if ( identical( class(res), "try-error") )  b2 <- b1
  m <- cbind(1, m1)
  m <- m / Rfast::rowsums(m)
  loglik <- sum( y * log(m), na.rm = TRUE )
  colnames(b2) <- paste("Y", 1:d, sep = "")
  rownames(b2) <- colnames(X)

  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- model.matrix( ~., data.frame(xnew) )
    if ( !con )  xnew <- xnew[, -1, drop = FALSE]
    mu <- cbind( 1, exp(xnew %*% b2) )
    est <- mu/Rfast::rowsums(mu)
    colnames(est) <- colnames(y)
  }

  list(iters = k, loglik = loglik, be = b2, est = est)
}



