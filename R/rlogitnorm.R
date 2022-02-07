rlogitnorm <- function(n, m, s, fast = FALSE){
  if (fast)  {
    x <- Rfast::Rnorm(n, m, s)
  } else   x <- rnorm(n, m, s)
  1 / (1 + exp(-x) )
}
