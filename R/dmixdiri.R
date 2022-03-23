dmixdiri <- function(x, a, prob, logged = TRUE) {
  g <- dim(a)[1]
  f <- 0
  for (i in 1:g)  f <- f + prob[i] * Compositional::ddiri(x, a[i, ], logged = FALSE)

  if ( logged )  f <- log(f)
  f
}
