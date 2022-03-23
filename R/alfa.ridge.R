alfa.ridge <- function(y, x, a, lambda, B = 1, xnew = NULL) {
  z <- Compositional::alfa(x, a, h = TRUE)$aff 
  mod <- Compositional::ridge.reg(y, z, lambda, B = B, xnew = xnew)
  mod 
}
