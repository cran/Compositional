el.test1 <- function(x, mu, R = 1, ncores = 1, graph = FALSE) {

  x <- as.matrix(x)
  mu <- as.vector(mu)

  runtime <- proc.time()
  res <- emplik::el.test(x, mu, maxit = 1000)
  runtime <- proc.time() - runtime

  if ( R <= 1 ) {
    res <- list(res = res, runtime = runtime)

  } else if (R > 1) {

    stat <- res$"-2LLR"
    tb <- numeric(R)
    n <- nrow(x)
    d <- ncol(x)
    m <- as.vector( Rfast::colmeans(x) )
    y <- x - rep( m - mu, rep(n, d) ) ## brings the data under the null hypothesis

    if (ncores == 1) {
      ## under the null hypothesis, i.e. mean vector equal to M
      runtime <- proc.time()
      for (i in 1:R) {
        b <- sample(1:n, n, replace = TRUE)
        tb[i] <- emplik::el.test(y[b, ], mu, maxit = 1000)$"-2LLR"
      }
      runtime <- proc.time() - runtime

    } else {
      runtime <- proc.time()
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl) ## make the cluster

      ww <- foreach( i = 1:R, .combine = cbind, .packages = "emplik",
        .export = "el.test" ) %dopar% {
         b <- sample(1:n, n, replace = TRUE)
         tb[i] <- emplik::el.test(y[b, ], mu, maxit = 1000)$"-2LLR"
     }

      stopCluster(cl) ## stop the cluster
      tb <- as.vector(ww)
      runtime <- proc.time() - runtime

    }

    pvalue <- ( sum(tb > stat) + 1 ) / (R + 1)
    res <- list(res = res, boot.pvalue = pvalue, runtime = runtime)

  }

  res

}




