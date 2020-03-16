totvar <- function(x, a = 0) {
  y <- Compositional::alfa(x, a, h = FALSE)$aff
  s <- Rfast::colVars(y)
  sum(s)
}
