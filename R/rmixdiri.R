rmixdiri <- function(n, a, prob) {
  p2 <- c(0, cumsum(prob))
  d <- dim(a)[2]  ## dimensionality of the data
  u <- runif(n)
  g <- dim(a)[1]  ## how many clusters are there
  ina <- as.numeric( cut(u, breaks = p2) )  ## the cluster of each observation
  ina <- sort(ina)
  nu <- as.vector( table(ina) )  ## frequency table of each cluster
  x <- NULL
  for (j in 1:g)  x <- rbind(x, Compositional::rdiri( nu[j], a[j, ] ) )
  list(id = ina, x = x)
}
