################################
#### Spatial sign covariance matrix
#### Tsagris Michail 2/2015
#### References: A Durre, D Vogel, DE Tyler (2014)
#### The spatial sign covariance matrix with unknown location
#### Journal of Multivariate Analysis, 130: 107--117.
#### http://arxiv.org/pdf/1307.5706v2.pdf
#### mtsagris@yahoo.gr
################################

sscov <- function(x, me = NULL, tol = 1e-09) {
  ## x contains the data

  x <- as.matrix(x)  ## makes sure x is a matrix
  n <- dim(x)[1]  ## sample size
  p <- dim(x)[2]

  if ( is.null(me) ) {
    me <- spat.med(x, tol)  ## spatial median of x
  } 

  y <- x - rep( me, rep(n, p) )
  rs <- sqrt ( Rfast::rowsums(y^2) )
  y <- y / rs  ## unit vectors
  crossprod( y ) / n  ## SSCM

}