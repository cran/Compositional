################################
#### Multivariate analysis of variance (MANOVA)
#### Tsagris Michail 8/2015
#### mtsagris@yahoo.gr
#### References: Johnson and Wichern (2007, 6th Edition)
#### Applied Multivariate Statistical Analysis p. 302-303.
###########
#### Todorov, Valentin and Filzmoser, Peter (2010)
#### Robust Statistic for the One-way MANOVA
#### Computational Statistics \& Data Analysis 54(1):37-48
################################

maov <- function(x, ina) {
  ## x is a matrix with the data
  ## ina is a grouping  variable indicating the groups

  x <- as.matrix(x)  ## makes sure x is a matrix
  ina <- as.numeric(ina)
  ni <- as.vector( table(ina) )  ## group sample sizes
  n <- nrow(x)  ## total sample size
  g <- max(ina)  ## number of groups
  p <- ncol(x)  ## dimensionality of the data
  s <- array( dim = c(p, p, g) )

  for (i in 1:g)  s[, , i] <- (ni[i] - 1) * cov(x[ina == i, ])  ## group covariances
  W <- apply(s, 1:2, sum)
  m <- aggregate(x, by = list(ina), mean)  ## group mean vectors
  m <- as.matrix(m[, -1])
  me <- colMeans(x)  ## total mean vector
  y <- sqrt(ni) * (m - rep(me, rep(g, p)) )
  B <- crossprod(y)
  A <- diag(p) + solve(W, B)
  lam <- 1/det(A)

  if (g == 2 ) {
    stat <- (n - p - 1 ) / p * (1 - lam)/lam
    pvalue <- pf( stat, p, n - p - 1, lower.tail = FALSE )
    note <- paste("F approximation has been used")

  } else if (g == 3) {
    stat <- (n - p - 2 )/p * (1 - sqrt(lam)) / sqrt(lam)
    pvalue <- pf( stat, 2 * p, 2 * (n - p - 2), lower.tail = FALSE )
    note <- paste("F approximation has been used")

  } else {
    stat <- -( n - 1 - (p + g)/2 ) * log(lam)
    pvalue <- pchisq( stat, p * (g - 1), lower.tail = FALSE )
    note <- paste("Chi-square approximation has been used")
  }

  result <- c(stat, pvalue)
  names(result) <- c('stat', 'p-value')
  list(note = note, result = result)

}
