################################
#### Normal mixture models for compositional data
#### Tsagris Michail 5/2015
#### mtsagris@yahoo.gr
#### References: Ryan P. Browne, Aisha ElSherbiny and
#### Paul D. McNicholas (2015)
#### R package mixture: Mixture Models for Clustering and Classification
################################
mix.compnorm <- function(x, g, model = NULL, type = "alr", veo = FALSE) {
  ## x is the compositional data
  ## g is the number of components to be used
  ## model is the type of model to be used
  ## type is either 'alr' or 'ilr'
  p <- dim(x)[2]  ## dimensionality of the data
  n <- dim(x)[1]  ## sample size

  if (type == "alr") {
    y <- log(x[, -1]/x[, 1])
  } else {
    y0 <- log(x)
    y1 <- y0 - Rfast::rowmeans( y0 )
    y <- tcrossprod( y1, helm(p) )
  }

  mod <- mixture::gpcm(y, G = g, mnames = model, start = 0, mmax = 100)
  param <- mod$gpar
  mu <- matrix(nrow = g, ncol = length( param[[ 1 ]]$mu) )
  su <- array( dim = c( length(param[[ 1 ]]$mu), length(param[[ 1 ]]$mu), g ) )

  for ( i in 1:g ) {
    mu[i, ] <- param[[ i ]]$mu  ## mean vector of each component
    su[, , i] <- param[[ i ]]$sigma  ## covariance of each component
  }
  prob <- param$pi  ## mixing probability of each component

  pij <- mod$z
  est <- Rfast::rowMaxs(pij)
  list(type = type, mu = mu, su = su, prob = prob, est = est)
}
