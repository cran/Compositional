################################
#### Flexible discriminant analysis for compositional data using the alpha-transformation
################################
alfa.fda <- function(xnew, x, ina, a) {
  y <- Compositional::alfa(x, a)$aff  ## apply the alpha-transformation
  xnew <- as.matrix(xnew)
  xnew <- matrix( xnew, ncol = dim(x)[2] )
  ynew <- Compositional::alfa(xnew, a)$aff
  mod <- mda::fda(ina ~ y)
  est <- predict(mod, ynew)
  list(mod = mod, est = est)
}
