lcreg.aov <- function(mod0, mod1) {
  sse0 <- sum( mod0$residuals^2 )
  sse1 <- sum( mod1$residuals^2 )
  p <- length(mod1$be)
  n <- length(mod1$residuals)
  stat <- (sse0 - sse1) / ( sse1 / (n - p) )
  pvalue <- pf(stat, 1, n - p, lower.tail = FALSE)
  res <- c(stat, pvalue)
  names(res) <- c( "F statistic", "p-value" )
  res
}
