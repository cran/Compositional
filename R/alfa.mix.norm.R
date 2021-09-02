################################
#### Normal mixture models for compositional data
#### Tsagris Michail 5/2015
#### mtsagris@yahoo.gr
#### References: Ryan P. Browne, Aisha ElSherbiny and
#### Paul D. McNicholas (2015)
#### R package mixture: Mixture Models for Clustering and Classification
################################
alfa.mix.norm <- function(x, g, a, model, veo = FALSE) {
  ## x is the compositional data
  ## g is the number of components to be used
  ## model is the type of model to be used
  ## type is either 'alr' or 'ilr'
  p <- dim(x)[2]  ## dimensionality of the data
  n <- dim(x)[1]  ## sample size

  z <- Compositional::alfa(x, a)
  y <- z$aff
  ja <- z$sa

  mod <- mixture::gpcm(y, G = g, mnames = model, start = 0, mmax = 100, veo = veo)
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
  list(mu = mu, su = su, prob = prob, est = est)
}
