########################
#### Multivariate kernel density estimation: tuning the bandwidth h
#### using the maximum likelihood cross validation
#### Tsagris Michail 2/2015
#### mtsagris@yahoo.gr
#### References: Arsalane Chouaib Guidoum (2015)
#### Kernel Estimator and Bandwidth Selection for Density
#### and its Derivatives. The kedd package
#### http://cran.r-project.org/web/packages/kedd/vignettes/kedd.pdf
#### M.P. Wand and M.C. Jones (1995)
#### Kernel smoothing, pages 91-92.
#### B.W. Silverman (1986)
#### Density estimation for statistics and data analysis, pages 76-78.
################################
mkde.tune <- function( x, low = 0.1, up = 3, s = cov(x) ) {
  ## x contains the multivariate data
  ## low and up are the limits within which the
  ## search is conducted
  n <- dim(x)[1]
  d <- dim(x)[2]  ## sample and dimensionality of x
  s <- s ## can put a robust covariance matrix here if you want
  eig <- eigen(s)
  lam <- eig$values  ## eigenvalues of the covariance matrix
  vec <- eig$vectors  ## eigenvectors of the covariance matrix
  B <- vec %*% ( t(vec) / sqrt(lam) )
  z <- x %*% B
  a2a <- Rfast::Dist( z, method = "euclidean", square = TRUE )
  a2a <- exp(-0.5 * a2a)
  ds <- 1 / prod(lam)^0.5

  tune <- function(h, n1) {
    a <- a2a^( 1 / h^2 )
    - d * log(h) + sum( log( rowSums( a ) - 1 ) ) / n1
  }

  low <- low     ;    up <- up
  bar <- optimize(tune, c(low, up), n1 = n - 1, maximum = TRUE)
  list( hopt = bar$maximum, maximum = bar$objective + log(ds) - d/2 * log(2 * pi) - log(n - 1) )
}
