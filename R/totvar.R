totvar <- function(x, a = 0) {
  y <- alfa(x, a, h = FALSE)$aff
  s <- Rfast::colVars(y)
  sum(s)
}
