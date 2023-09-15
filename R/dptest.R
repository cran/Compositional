dptest <- function(x1, x2, B = 100) {
  d <- dim(x1)[2]
  u <- Directional::riag(B, numeric(d))
  x1 <- sqrt(x1)  ;   x2 <- sqrt(x2)
  z1 <- tcrossprod(x1, u) 
  z2 <- tcrossprod(x2, u) 
  p <- numeric(B)
  for (i in 1:B)  
  p[i] <- ks.test(z1[, i], z2[, i])$p.value
  p <- sort(p)
  i <- 1:B
  pval <- min(B / i * p )
  list(pvalues = p, pvalue = pval)
}

