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

mkde <- function(x, h, thumb = "none") {
  x <- as.matrix(x)  ## makes sure x is a matrix
  ## h is the h you want, which is either a vector or a single number
  ## thumb can be either "none" so the specified h is used, or
  ## "scott", or "silverman"
  n <- nrow(x)
  d <- ncol(x)  ## sample and dimensionality of x

  if ( thumb == "silverman" ) {
    s <- Rfast::colVars(x, std = TRUE) 
    h <- ( 4/(d + 2) )^( 1/(d + 4) ) * s * n^( -1/(d + 4) )

  } else  if ( thumb == "scott" ) {
    s <- Rfast::colVars(x, std = TRUE) 
    h <- s * n^( -1/(d + 4) )

  } else if ( thumb == "estim" ) {
    h <- mkde.tune(x)$hopt

  } else  h <- h

  if ( length(h) == 1 ) {
    h <- diag( 1 / h, d )
  } else h <- diag( 1 / h)

  con <- prod( diag( h ) )
  y <- x %*% h
  a1 <- fields::rdist(y, compact = FALSE)
  f <- 1/(2 * pi)^(d/2) * con * as.vector( Rfast::rowmeans( exp(-0.5 * a1^2 ) ) )
  f

}
