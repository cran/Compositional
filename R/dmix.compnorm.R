dmix.compnorm <- function(x, mu, sigma, prob, type = "alr", logged = TRUE) {

  if ( is.null(dim(x)[1]) )  x <- matrix(x, nrow = 1)
  g <- dim(mu)[1]

  if (type == "alr") {
    y <- Compositional::alr(x) # additive log-ratio transformation
  } else if (type == "ilr") {
    y <- Compositional::alfa(x, 0)$aff
  } else  y <- Compositional::pivot(x)

  f <- 0
  for (i in 1:g)  f <- f + prob[i] * Rfast::dmvnorm(y, mu[i, ], sigma[, , i])

  if ( logged )  f <- log(f)
  f

}
