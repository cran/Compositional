maov <- function(x, ina) {
  ina <- as.numeric(ina)
  ni <- tabulate(ina)  ## group sample sizes
  n <- dim(x)[1]  ## total sample size
  g <- max(ina)  ## number of groups
  p <- dim(x)[2]  ## dimensionality of the data
  m <- Rfast2::colGroup(x, as.integer(ina) ) / ni
  me <- Rfast::colmeans(x)  ## total mean vector
  y <- sqrt(ni) * (m - rep(me, rep(g, p)) )
  B <- crossprod(y)
  Tot <- Rfast::cova(x) * (n - 1)
  lam <- det(Tot - B) / det(Tot)

  if ( g == 2 ) {
    stat <- (n - p - 1 ) / p * (1 - lam)/lam
    pvalue <- pf( stat, p, n - p - 1, lower.tail = FALSE )
    note <- paste("F approximation has been used")
  } else if ( g == 3 ) {
    stat <- (n - p - 2 )/p * ( 1 - sqrt(lam) ) / sqrt(lam)
    pvalue <- pf( stat, 2 * p, 2 * (n - p - 2), lower.tail = FALSE )
    note <- paste("F approximation has been used")
  } else {
    stat <-  -( n - 1 - (p + g)/2 ) * log(lam)
    pvalue <- pchisq( stat, p * (g - 1), lower.tail = FALSE )
    note <- paste("Chi-square approximation has been used")
  }

  result <- c(stat, pvalue)
  names(result) <- c('stat', 'p-value')
  list(note = note, result = result)
}
