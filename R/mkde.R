########################
#### Multivariate kernel density estimation
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
mkde <- function(x, h = NULL, thumb = "silverman") {
  ## h is the h you want, which is either a vector or a single number
  ## thumb can be either "none" so the specified h is used, or
  ## "scott", or "silverman"
  n <- dim(x)[1]
  d <- dim(x)[2]  ## sample and dimensionality of x

  if ( is.null(h) ) {
    if ( thumb == "silverman" ) {
      s <- Rfast::colVars(x, std = TRUE)
      h <- ( 4/(d + 2) )^( 1/(d + 4) ) * s * n^( -1/(d + 4) )
    } else  if ( thumb == "scott" ) {
      s <- Rfast::colVars(x, std = TRUE)
      h <- s * n^( -1/(d + 4) )
    } else if ( thumb == "estim" ) {
      h <- Compositional::mkde.tune(x)$hopt
    }
  }

  if ( length(h) == 1 ) {
    h <- diag( 1 / h, d )
  } else h <- diag( 1 / h)

  con <- prod( diag( h ) )
  y <- x %*%  h
  a1 <- Rfast::Dist(y, method = "euclidean", square = TRUE)
  (0.5 / pi)^(d/2) * con * ( Rfast::rowsums( exp( - 0.5 * a1 ) ) - 1 ) / (n - 1)
}
