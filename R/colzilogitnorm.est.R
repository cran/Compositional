colzilogitnorm.est <- function(x) {
  p <- dim(x)[2]
  param <- matrix(nrow = p, ncol = 4)
  colnames(param) <- c("prop1", "mean", "unbiased variance", "loglik")
  for (i in 1:p) {
    mod <- Compositional::zilogitnorm.est(x[, i])
    par[i, ] <- c(mod$param, mod$loglik)
  }
  param
}
